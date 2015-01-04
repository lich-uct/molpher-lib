/*
 Copyright (c) 2012 Petr Škoda

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <algorithm>
// remove
#include <iostream>
#include <math.h>
#include <limits>
#include <fstream>
#include <sstream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"

#include "PcaReducer.h"

//#define LOG_PCA_DATA

/**
 * Return false if given value is nan of inf.
 */
bool inline ValidDouble(const double& value)
{
    return !(
        (value == std::numeric_limits<double>::infinity())||
        (value != value) );
}

/**
 * Search for original file name by adding numbers.
 */
int FileNumber(const std::string& prefix, const std::string& suffix)
{
    int number = 1;
    bool found = false;
    
    while(true)
    {
        std::stringstream ss;
        ss << prefix << number << suffix;
        
        std::string fileName = ss.str();
        // try to open file
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::out | std::ios_base::in);
        if (file.is_open())
        { // file exist
            file.close();
            ++number;
        }
        else
        {
            break;
        }
    } 
    return number;
}

/**
 * Return square of maximum value of error which is feasible when
 * caltulating eigen vectors. 
 * @return double Square of maximum error.
 */
double inline MaxError() { return 0.00001 * 0.00001; }

/**
 * Multiply square matrix with vector. Matrix
 * is multiplied by vector from right.
 * Vector and result must not be the same.
 * Values in matrix  must be stored as rows.
 * @param[in] const double* vector Vector.
 * @param[in] double* result Result of operation.
 * @param[out] const double* matrix Square matrix.
 * @param[in] size_t dataDimension Size of vector or matrix.
 */
void inline Multiply(
        const double* vector, 
        double* result,
        const double* matrix, 
        size_t dataDimension) {
    // store output vector in result
    for (int i = 0; i < dataDimension; ++i) {
        // calculate start line index
        int lineIndex = i * dataDimension;
        // set initial value to zero
        result[i] = 0;
        for (int j = 0; j < dataDimension; ++j, ++lineIndex) {                
            result[i] += vector[j] * matrix[lineIndex];            
            assert(ValidDouble(result[i]));
        }
    }
}

/**
 * Calculate error between two vector.
 * Error is calculated as powered Euclidean distance.
 * @param[in] const double* oldVector Vector.
 * @param[in] const double* newVector Vector.
 * @param[in] size_t dataDimension Size of vectors.
 * @return double
 */
double inline CalculateError(const double* oldVector, const double* newVector, size_t dataDimension) {
    // according to formule : sum(pow(a[i]-b[i], 2)
    double result = 0;
    for (int i = 0; i < dataDimension; ++i){
        result += std::pow(oldVector[i] - newVector[i], 2);
    }
    assert(ValidDouble(result));
    return result;
}

/**
 * Calculates scalar multiplies of vectors
 * @param[in] const double* left Vector.
 * @param[in] const double* right Vector.
 * @param[in] size_t dataDimension Size of vectors.
 * @return double Scalar of two vectors.
 */
double inline MultiplyScalar(const double* left, const double* right, size_t dataDimension) {
	// according to formule : s = sum( left.xi * right.xi] )
    double result = 0;
    for (int i = 0; i < dataDimension; ++i){
        result += left[i] * right[i];
    }
    assert(ValidDouble(result));
    return result;
}

/**
 * Normalize given vector.
 * @param[in,out] double* Vector to normalize.
 * @param[in] size_t dataDimension Vector size.
 */
void inline Normalize(double* vector, size_t dataDimension) {
    // calculate size of original vector
    double vectorSize = 0;
    for (int i = 0; i < dataDimension; ++i) {
        vectorSize += vector[i] * vector[i];
    }
    
    if (vectorSize == 0) {
        // input is zero vector .. 
        assert(false);
        return;
    }

    vectorSize = std::sqrt(vectorSize);
    // normalize
    for (int i = 0; i < dataDimension; ++i) 
    {
        vector[i] /= vectorSize;
        assert(ValidDouble(vector[i]));
    }
}

/**
 * Makes two vectors ortogonal. First given vector 
 * is fixed, second vector is changed to be ortogonal to first.
 * Result vector is not normalized.
 * @param[in] const double* first Normalized vector.
 * @param[in] const double* second Vector to ortogonalized.
 * @param[out] double* result Ortogonalized second vector.
 * @param[in] dataDimension Vectors size.
 */
void inline Orthogonalized(
        const double* first,
        const double* second,
        double* result,
        size_t dataDimension) {
    /*
     * use fallowing alhorihm:
     * alpha = <first,second>
     * result = second - alpha * first;
     */
    double alpha = MultiplyScalar(first, second, dataDimension);
    for (int i = 0; i < dataDimension; ++i) {
        result[i] = second[i] - (alpha * first[i]);
        assert(ValidDouble(result[i]));
    }
}

PcaReducer::PcaReducer() {
}

PcaReducer::~PcaReducer() {
}

void PcaReducer::Reduce(
        MolPtrVector &mols,
        FingerprintSelector fingerprintSelector,
        SimCoeffSelector simCoeffSelector,
        tbb::task_group_context &tbbCtx) {
    
    // check if we have some input data, 
    // alse prevent division by zero when calculating coordinates mean
    if (mols.size() == 0 || mols.size() == 1) {
        return;
    }
#ifdef LOG_PCA_DATA    
int fileNumber = FileNumber("m", ".txt");
#endif

    MeasureStage measureStage;
    SimCoefCalculator calc(simCoeffSelector, fingerprintSelector);
    // sotre number of objects == mols
    size_t objectsCount = mols.size();
    // caltulate fingerprints for all molecules
    std::vector<Fingerprint *> fingerprints;
    fingerprints.resize(objectsCount, NULL);
    CalculateFingerprints calculateFingerprints(calc, mols, fingerprints);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mols.size()),
            calculateFingerprints, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateFingerprints");
    }    
    
    // we calculate number of coordinates, and allocate them in a single new    
    // each FP will be split into units
    int coordinatesDimension = fingerprints[0]->getNumBits() / 32;
    // create space for coordinates, store data in single uber-array    
    size_t coordinatesSize = objectsCount * coordinatesDimension;
    double *coordinates = new double[coordinatesSize];
    // now we need transform FP into coordinates -> we use SUM method (parallel)
    CalculateCoordinatesSum calculateCoordinatesSum(fingerprints, coordinates, coordinatesDimension);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, objectsCount),
            calculateCoordinatesSum, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateCoordinates(Sum)");
    }    
    // we dont need FP no more
    for (size_t i = 0; i < mols.size(); ++i) {
        delete fingerprints[i];
    }
    fingerprints.clear();
    // center input coordinates
    
    double *meanCoordinates = new double[coordinatesDimension];
    // calculate mean
    if (!Cancelled(tbbCtx)) {
        // set to zero
        for (int j = 0; j < coordinatesDimension; ++j) {
            meanCoordinates[j] = 0;
        }        
        // we go throug data
        for (int i = 0; i < objectsCount; ++i) {
            // determine start coordinates
            size_t coordinateIndex = i * coordinatesDimension;
            for (size_t j = 0; j < coordinatesDimension; ++j, ++coordinateIndex) {
                // add to mean 
                meanCoordinates[j] += coordinates[coordinateIndex];
            }
        }
        for (int j = 0; j < coordinatesDimension; ++j) {
            meanCoordinates[j] /= (double)objectsCount;
        }
        measureStage.ReportAndReset("MeanCalculated");
    }    
    
    // center data (parallel)
    CenterCoordinates centerCoordinates(coordinates, meanCoordinates, coordinatesDimension);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, objectsCount),
            centerCoordinates, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CenterCoordinates");
    }    
    // delete mean
    delete[] meanCoordinates;
    meanCoordinates = 0;
    
#ifdef LOG_PCA_DATA     
{
    std::stringstream ss;
    ss << "coord" << fileNumber << ".txt";
    std::ofstream outfile(ss.str().c_str(), std::ios::out | std::ios::app);    
    // outfile << "covariance matrix" << std::endl;
    for (int o = 0; o < objectsCount; ++o) {    
        outfile << coordinates[o * coordinatesDimension];
        for (int j = 1; j < coordinatesDimension; ++j) {
            outfile << "," << coordinates[j + (o * coordinatesDimension)];
        }
        outfile << std::endl;
    }
    
    outfile.close();
}    
measureStage.ReportAndReset("LogOutput");    
#endif

    // compute covariance matrix, we can use the fact that it is symetric (paralel)
    double *covarianceMatrix = new double[coordinatesDimension * coordinatesDimension];
    CalculateCovarianceMatrix calCovMatrix(
        coordinates, covarianceMatrix, coordinatesDimension, objectsCount);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range2d<size_t, size_t>(0, coordinatesDimension, 0, coordinatesDimension),
            calCovMatrix, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateCovarianceMatrix");
    }    

    // now we need 2 most significant eigen vectors
    double* eigenFirst = new double[coordinatesDimension];
    double* eigenSecond = new double[coordinatesDimension];
    double* tempVector = new double[coordinatesDimension];
    // set vector to ones
    for (int i = 0; i < coordinatesDimension; ++i) {
        eigenFirst[i] = eigenSecond[i] = 1;
    }
    
    // calculate eigens ..     
    if (!Cancelled(tbbCtx)) {
        int iter = 0;
        // start with first eigen vector
        do {
            ++iter;
                     
            // multiply vector with matrix and store result into tempVector
            Multiply(eigenFirst, tempVector, covarianceMatrix, coordinatesDimension);
            // move data from tempVector into eigenFirst
            std::swap(eigenFirst, tempVector);
            // normalize before caltulating error
            Normalize(eigenFirst, coordinatesDimension); 
            Normalize(tempVector, coordinatesDimension); 
        } 
        while ( CalculateError(eigenFirst, tempVector, coordinatesDimension) > MaxError() );
        measureStage.ReportAndReset("CalculateFirstEigenVector");    
        // eigenFirst is already normalized
    }     
    
    if (!Cancelled(tbbCtx)) {
        // now we need second eigen vector, eigen vectors are ortogonal    
        do {
            // multiply vector with matrix and store result into tempVector
            Multiply(eigenSecond, tempVector, covarianceMatrix, coordinatesDimension);
            // make eigenSecond ortogonal to eigenFirst
            Orthogonalized(eigenFirst, tempVector, eigenSecond, coordinatesDimension);
            // normalize before caltulating error
            Normalize(eigenSecond, coordinatesDimension); 
            Normalize(tempVector, coordinatesDimension);             
        } 
        while ( CalculateError(eigenSecond, tempVector, coordinatesDimension) > MaxError() );
        measureStage.ReportAndReset("CalculateSecondEigenVector");
        // eigenSecond is already normalized
    }
#ifdef LOG_PCA_DATA     
{
    std::stringstream ss;
    ss << "m" << fileNumber << ".txt";
    std::ofstream outfile(ss.str().c_str(), std::ios::out | std::ios::app);    
    // outfile << "covariance matrix" << std::endl;
    for (int y = 0; y < coordinatesDimension; ++y) {
    outfile << covarianceMatrix[y * coordinatesDimension];
    for (int x = 1; x < coordinatesDimension; ++x) {
        outfile << "," << covarianceMatrix[x+(y * coordinatesDimension)];
    }
    outfile << std::endl;
    }
    outfile.close();
}

{
    std::stringstream ss;
    ss << "e1" << fileNumber << ".txt";
    std::ofstream outfile(ss.str().c_str(), std::ios::out | std::ios::app);
    // outfile << "eigen first" << std::endl;
    outfile << eigenFirst[0];
    for (int i = 1; i < coordinatesDimension; ++i) {
        outfile << "," << eigenFirst[i];
    }
    outfile << std::endl;
    outfile.close();
}

{
    std::stringstream ss;
    ss << "e2" << fileNumber << ".txt";
    std::ofstream outfile(ss.str().c_str(), std::ios::out | std::ios::app);
    // outfile << "eigen second" << std::endl;
    outfile << eigenSecond[0];
    for (int i = 1; i < coordinatesDimension; ++i) {
        outfile << "," << eigenSecond[i];
    }
    outfile << std::endl;
    outfile.close();
}   

measureStage.ReportAndReset("LogOutput");
#endif
    // delete covariance matrix
    delete[] covarianceMatrix;
    covarianceMatrix = 0;    

    // on the end we transform data (parallel)
    CalculateCoordinates calculateCoordinates(coordinates, 
            eigenFirst, eigenSecond, mols, coordinatesDimension);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, objectsCount),
            calculateCoordinates, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateCoordinates");
    }

    // delete coordinates
    delete[] coordinates;
    coordinates = 0;
    // delete eigens
    delete[] eigenFirst;
    eigenFirst = 0;
    delete[] eigenSecond;
    eigenSecond = 0;
    delete[] tempVector;
    tempVector = 0;    
}

bool PcaReducer::Cancelled(tbb::task_group_context &ctx)
{
    return ctx.is_group_execution_cancelled();
}

PcaReducer::CalculateFingerprints::CalculateFingerprints(SimCoefCalculator &calc,
            MolPtrVector &mols, std::vector<Fingerprint *> &fingerprints)
        : mCalc(calc), mMols(mols), mFingerprints(fingerprints)
{ }

void PcaReducer::CalculateFingerprints::operator()(
    const tbb::blocked_range<size_t> &r) const
{
    for (size_t i = r.begin(); i != r.end(); ++i) {
        RDKit::RWMol *mol = NULL;
        try {
            mol = RDKit::SmilesToMol(mMols[i]->smile);
            if (mol) {
                RDKit::MolOps::Kekulize(*mol);
                mFingerprints[i] = mCalc.GetFingerprint(mol);
                delete mol;
            } else {
                throw ValueErrorException("");
            }
        } catch (const ValueErrorException &exc) {
            mFingerprints[i] = NULL;
            delete mol;
        }
    }
}

PcaReducer::MeasureStage::MeasureStage()
{
    mTimestamp = std::clock();
}

void PcaReducer::MeasureStage::ReportAndReset(const std::string &stage)
{
    clock_t current = std::clock();
#if PCAREDUCER_REPORTING == 1
    std::ostringstream stream;
    stream << "PcaReducer: " << stage <<
        " consumed " << current - mTimestamp << " msec.";
    SynchCout(stream.str());
#endif
    mTimestamp = current;
}

PcaReducer::CalculateCoordinatesSum::CalculateCoordinatesSum(
        const std::vector<Fingerprint *>& fingerprints, double *coordinates,  size_t outDimension)
        : mFingerprints(fingerprints), mCoordinates(coordinates), mOutDimension(outDimension)
{ }

void PcaReducer::CalculateCoordinatesSum::operator()(
        const tbb::blocked_range<size_t> &param) const 
{
    for (size_t f = param.begin(); f != param.end(); ++f) 
    {
       // we calculate coordinates for i-th FP
       Fingerprint* fp = mFingerprints[f];        
       // calculate size of single interpret unit
       // also we want round result up 
       int unitSizeMax = (int)
         ( ( (double)fp->getNumBits() / mOutDimension) + 0.5 );
       // check that unit size is greater than zero
       assert(unitSizeMax > 0);        
       // for each unit in FP
       size_t coordStart = f * mOutDimension; // save start position for coordinates
       size_t coordMax = (f + 1) * mOutDimension;
       // used to determine size of interpret unit
       size_t unitSize = 0;
       // bite position
       int bitePos = 0;
       int biteMax = fp->getNumBits();
       for (size_t c = coordStart; c < coordMax; ++c) 
       {
           mCoordinates[c] = 0;
           for (size_t i = 0; i < unitSizeMax && bitePos < biteMax; ++i) 
           {
               // store bite value
               mCoordinates[c] += (*fp->dp_bits)[bitePos];
               // increase bite
               ++bitePos;
           }
       }       
    }
}

PcaReducer::CenterCoordinates::CenterCoordinates(
        double *coordinates, const double *mean, size_t dimension)
        : mCoordinates(coordinates), mMean(mean), mDimension(dimension)
{}

void PcaReducer::CenterCoordinates::operator()(
        const tbb::blocked_range<size_t> &param) const
{
    // we substract mean from coordinates
    for (size_t f = param.begin(); f != param.end(); ++f) {
        size_t coordinatesStart = mDimension * f;
        for (size_t i = 0; i < mDimension; ++i, ++coordinatesStart) {
            mCoordinates[coordinatesStart] -= mMean[i];
        }
    }
}

PcaReducer::CalculateCovarianceMatrix::CalculateCovarianceMatrix(
        const double *coordinates, double *matrix, size_t dimension, size_t objectCount)
        : mCoordinates(coordinates), mMatrix(matrix), mDimension(dimension), 
          mObjectCount(objectCount)
{ }

void PcaReducer::CalculateCovarianceMatrix::operator()(
        const tbb::blocked_range2d<size_t, size_t> &param) const
{
//    int coordinatesSize = mObjectCount * mDimension;
    for (size_t r = param.rows().begin(); r != param.rows().end(); ++r) {
        for (size_t c = param.cols().begin(); c != param.cols().end(); ++c) {
            // we should caldulate for i,j 
            // use symetricity, calculate only for i >= j
            if (r < c) break;
            // otherwise start computing ...
            double result = 0;
            // i jump on begin of data
            size_t pos = 0;
            for (size_t i = 0; i < mObjectCount; ++i, pos += mDimension) {
                // i + r and i + c point to r,c-th value of i-th object
                result += mCoordinates[pos + r] * mCoordinates[pos + c];
            }
            // divide result with mObjectCount
            result /= (mObjectCount - 1);
            // use matrix symetricity and store result
            mMatrix[r + (c * mDimension)] = 
                    mMatrix[c + (r * mDimension)] = 
                        result;
        }   
    }
}

PcaReducer::CalculateCoordinates::CalculateCoordinates(
                const double *coordinates, const double *eigenFirst, const double* eigenSecond, 
                MolPtrVector& mols, size_t dimension)
        : mCoordinates(coordinates), mEigenFirst(eigenFirst), 
          mEigenSecond(eigenSecond), mMols(mols), mDimension(dimension)
 { }
void PcaReducer::CalculateCoordinates::operator()(
        const tbb::blocked_range<size_t> &param) const
{
    // calculate coordinates, multiply original coordinates with eigens
    // to gain result coordinates
    for (size_t f = param.begin(); f != param.end(); ++f) {
        // start by setting posX, ans pos Y to zero
        mMols[f]->posX = 0;
        mMols[f]->posY = 0;
        // we calculate start index, in mCoordination for object with index f
        int coordPos = f * mDimension;
        // multiply vectors
        for (int i = 0; i < mDimension; ++i, ++coordPos) {
            mMols[f]->posX += mCoordinates[coordPos] * mEigenFirst[i];
            mMols[f]->posY += mCoordinates[coordPos] * mEigenSecond[i];
        }
    }
}

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

#pragma once

#include <string>
#include <vector>
#include <ctime>

#include <tbb/task.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_reduce.h>

#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "DimensionReducer.h"
#include "chem/simCoefStrategy/SimCoefStrategy.h"
#include "chem/SimCoefCalculator.hpp"

#ifndef PCAREDUCER_REPORTING
#define PCAREDUCER_REPORTING 1
#endif

class PcaReducer : public DimensionReducer
{
public:
    PcaReducer();
    ~PcaReducer();
    virtual void Reduce(
        MolPtrVector& mols,
        FingerprintSelector fingerprintSelector,
        SimCoeffSelector simCoeffSelector,
        tbb::task_group_context& tbbCtx);
protected:   
    bool Cancelled(tbb::task_group_context &ctx);    
protected:
    /**
     * Class is used in cooperation with tbb. Class compute
     * fingerprints according to given SimCoefCalculator calculator.
     */
    class CalculateFingerprints
    {
    public:
        /**
         * Base ctor.
         * @param SimCoefCalculator calc
         * @param MolPtrVector& mols
         * @param std::vector<Fingerprint *> fingerprints Fingerprints storage.
         */
        CalculateFingerprints(SimCoefCalculator &calc,
            MolPtrVector& mols, std::vector<Fingerprint *>& fingerprints);
        /**
         * Operator for tbb.
         */
        void operator()(const tbb::blocked_range<size_t>& r) const;
    private:
        /**
         * Molpher class used to genererate fingerprints.
         */
        SimCoefCalculator& mCalc;
        /**
         * Molecule storage.
         */
        MolPtrVector& mMols;
        /**
         * Fingerprints storage.
         */
        std::vector<Fingerprint *>& mFingerprints;
    };
    /**
     * Class is used to measure and report time, that
     * has been spend in subparts of Reduce method.
     */
    class MeasureStage
    {
    public:
        /**
         * Base ctor, set inner timer to current time.
         */
        MeasureStage();
        /**
         * Report string with time from timer and reset timer to currnt timer.
         * @param const std::string& stage String to report.
         */
        void ReportAndReset(const std::string& stage);
    private:
        /**
         * Store time of last report.
         */
        clock_t mTimestamp;
    };
    /**
     * Class transform fingerprints into coordinates. Class
     * is desing to be used with tbb.
     */
    class CalculateCoordinatesSum
    {
    public:
        /**
         * Base ctor.
         * @param std::vector<Fingerprint *>& fingerprints Fingeprtins.
         * @param coordinates double* Array with coordinates of objects.
         * @param outDimension size_t Store dimension of output coordinates.
         */
        CalculateCoordinatesSum(const std::vector<Fingerprint *>& fingerprints, 
                double* coordinates, 
                size_t outDimension);
        /**
         * Operator for tbb.
         */        
        void operator()(const tbb::blocked_range<size_t>& param) const;
    private:
        /**
         * Store fingerprints.
         */
        const std::vector<Fingerprint *>& mFingerprints;
        /**
         * Pointer to output coordinates storage.
         */
        double* mCoordinates;
        /**
         * Coordinates dimension.
         */
        size_t mOutDimension;
    };
    /**
     * Center coordinates = set. coordinates mean to zero.
     * Class is desing to be used with tbb.
     */
    class CenterCoordinates
    {
    public:
        /**
         * Base ctor.
         * @param double* coordinates Array with coordinates of objects.
         * @param const double* Mean used to centred data.
         * @param size_t dimension Data dimension (coordinates per object)
         */
        CenterCoordinates(double *coordinates, const double *mean, size_t dimension);
        /**
         * Operator for tbb.
         */          
        void operator()(const tbb::blocked_range<size_t>& param) const;
    private:
        /**
         * Coordinates storage.
         */
        double* mCoordinates;
        /**
         * Data mean.
         */
        const double* mMean;
        /**
         * Coordinates dimension.
         */
        size_t mDimension;
    };
    /**
     * Class calculate a covarinace matrix for given coordinates.
     * Class is desing to be used with tbb.
     */
    class CalculateCovarianceMatrix
    {
    public:
        /**
         * Base ctor.
         * @param const double* coordinates Array with coordinates of objects.
         * @param double* matrix Covariance matrix.
         * @param size_t dimension Object dimension (coordinates per object)
         * @param size_t objectCount Number of objects.
         */
        CalculateCovarianceMatrix(
                const double* coordinates, double* matrix, 
                size_t dimension, size_t objectCount);
        /**
         * Operator for tbb.
         */         
        void operator()(const tbb::blocked_range2d<size_t, size_t>& param) const;
    private:
        /**
         * Coordinates storage.
         */
        const double* mCoordinates;
        /**
         * Covariance matrix.
         */
        double* mMatrix;
        /**
         * Data dimension.
         */
        size_t mDimension;
        /**
         * Data count.
         */
        size_t mObjectCount;
    };
    /**
     * Transform coordinates with given vectors(eigenVectors) and
     * store result into molecule storage(mols). Class
     * is desinged to be used with tbb.
     */
    class CalculateCoordinates
    {
    public:
        /**
         * Basector.
         * @param double* coordinates  Data dimension (coordinates per object)
         * @param const double* eigenFirst Largets eigen vector.
         * @param const double* eigenSecond Second largest eigen vector.
         * @param mols MolPtrVector& Acces to molecule storage, used to store positions.
         * @param dimension  Data dimension (coordinates per object)
         */
        CalculateCoordinates(
                const double* coordinates, 
                const double* eigenFirst, const double* eigenSecond, 
                MolPtrVector& mols, size_t dimension);
        /**
         * Operator for tbb.
         */            
        void operator()(const tbb::blocked_range<size_t>& param) const;
    private:
        /**
         * Coordinates storage.
         */
        const double* mCoordinates;
        /**
         * First eigen vector.
         */
        const double* mEigenFirst;
        /**
         * Second eigen vector.
         */
        const double* mEigenSecond;
        /**
         * Molecule storage.
         */
        MolPtrVector& mMols;
        /**
         * Coordinates dimension.
         */
        size_t mDimension;        
    };    
};
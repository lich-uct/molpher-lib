/*
 Copyright (c) 2012 Petr Koupy

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

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"
#include "KamadaKawaiReducer.h"

KamadaKawaiReducer::KamadaKawaiReducer()
{
}

KamadaKawaiReducer::~KamadaKawaiReducer()
{
}

bool KamadaKawaiReducer::Cancelled(tbb::task_group_context &ctx)
{
    return ctx.is_group_execution_cancelled();
}

KamadaKawaiReducer::MeasureStage::MeasureStage()
{
    mTimestamp = std::clock();
}

void KamadaKawaiReducer::MeasureStage::ReportAndReset(const std::string &stage)
{
    clock_t current = std::clock();
#if KAMADAKAWAIREDUCER_REPORTING == 1
    std::ostringstream stream;
    stream << "KamadaKawaiReducer: " << stage <<
        " consumed " << current - mTimestamp << " msec.";
    SynchCout(stream.str());
#endif
    mTimestamp = current;
}

void KamadaKawaiReducer::RandomizeCoordinates::operator()(
    const tbb::blocked_range<MolPtrVector::iterator> &mols) const
{
    MolPtrVector::iterator it;
    for (it = mols.begin(); it != mols.end(); it++) {
        (*it)->posX = SynchRand::GetRandomNumber(0, 800);
        (*it)->posY = SynchRand::GetRandomNumber(0, 800);
    }
}

KamadaKawaiReducer::CalculateFingerprints::CalculateFingerprints(
    SimCoefCalculator &calc, MolPtrVector &mols,
    std::vector<Fingerprint *> &fingerprints
    ) :
    mCalc(calc),
    mMols(mols),
    mFingerprints(fingerprints)
{
    assert(mMols.size() == mFingerprints.size());
}

void KamadaKawaiReducer::CalculateFingerprints::operator()(
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

KamadaKawaiReducer::CalculateDistances::CalculateDistances(
    SimCoefCalculator &calc, MolPtrVector &mols,
    std::vector<Fingerprint *> &fingerprints, double **dist
    ) :
    mCalc(calc),
    mMols(mols),
    mFingerprints(fingerprints),
    mDist(dist)
{
    assert(mDist);
}

void KamadaKawaiReducer::CalculateDistances::operator()(
    const tbb::blocked_range2d<size_t> &r) const
{
    for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
        for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
            double dist = 1.0;
            if (mFingerprints[i] && mFingerprints[j]) {
                double coeff = mCalc.GetSimCoef(mFingerprints[i], mFingerprints[j]);
                bool connected = true;
                    //(mMols[i]->parentSmile == mMols[j]->smile) ||
                    //(mMols[j]->parentSmile == mMols[i]->smile);
                if (connected) {
                    dist = mCalc.ConvertToDistance(coeff);
                } else {
                    dist = 2 * mCalc.ConvertToDistance(coeff);
                }
            }
            mDist[i][j] = dist * 100;
        }
    }
}

KamadaKawaiReducer::IterateFloydWarshall::IterateFloydWarshall(
    size_t iter, double **last, double **current
    ) :
    mIter(iter),
    mLast(last),
    mCurrent(current)
{
    assert(mLast);
    assert(mCurrent);
}

void KamadaKawaiReducer::IterateFloydWarshall::operator()(
    const tbb::blocked_range2d<size_t> &r) const
{
    for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
        for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
            double oldDist = mLast[i][j];
            double newDist = mLast[i][mIter] + mLast[mIter][j];
            mCurrent[i][j] = newDist < oldDist ? newDist : oldDist;
        }
    }
}

KamadaKawaiReducer::CalculateStiffness::CalculateStiffness(
    double stiffnessFactor, double **dist, double **stiff
    ) :
    mFactor(stiffnessFactor),
    mDist(dist),
    mStiff(stiff)
{
    assert(mDist);
    assert(mStiff);
}

void KamadaKawaiReducer::CalculateStiffness::operator()(
    const tbb::blocked_range2d<size_t> &r) const
{
    for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
        for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
            double d = mDist[i][j];
            if (d > 0) {
                mStiff[i][j] = mFactor / (d * d);
            } else {
                mStiff[i][j] = 0;
            }
        }
    }
}

KamadaKawaiReducer::FindMaxGradient::FindMaxGradient(
    double separationFactor, MolPtrVector &mols, double **dist, double **stiff
    ) :
    mFactor(separationFactor),
    mMaxGradient(0),
    mMaxGradientIdx(0),
    mMols(mols),
    mDist(dist),
    mStiff(stiff)
{
    assert(mDist);
    assert(mStiff);
}

KamadaKawaiReducer::FindMaxGradient::FindMaxGradient(
    FindMaxGradient &toSplit, tbb::split
    ) :
    mFactor(toSplit.mFactor),
    mMaxGradient(0),
    mMaxGradientIdx(0),
    mMols(toSplit.mMols),
    mDist(toSplit.mDist),
    mStiff(toSplit.mStiff)
{
}

void KamadaKawaiReducer::FindMaxGradient::operator()(
    const tbb::blocked_range<size_t> &r)
{
    for (size_t i = r.begin(); i != r.end(); ++i) {
        double dEdx = 0; // dE/dx
        double dEdy = 0; // dE/dy

        for (int j = 0; j < mMols.size(); ++j) {
            double dx = mMols[i]->posX - mMols[j]->posX;
            double dy = mMols[i]->posY - mMols[j]->posY;
            double euclidDistance = sqrt(dx * dx + dy * dy);

            if (euclidDistance > 0) {
                double stiffness = mStiff[i][j];
                double distance = mDist[i][j];

                dEdx += stiffness * (dx - (mFactor * distance * dx) / euclidDistance);
                dEdy += stiffness * (dy - (mFactor * distance * dy) / euclidDistance);
            }
        }

        double currentGradient = sqrt(dEdx * dEdx + dEdy * dEdy);
        if (currentGradient > mMaxGradient) {
            mMaxGradient = currentGradient;
            mMaxGradientIdx = i;
        }
    }
}

void KamadaKawaiReducer::FindMaxGradient::join(FindMaxGradient &toJoin)
{
    if (toJoin.mMaxGradient > mMaxGradient) {
        mMaxGradient = toJoin.mMaxGradient;
        mMaxGradientIdx = toJoin.mMaxGradientIdx;
    }
}

double KamadaKawaiReducer::FindMaxGradient::GetMaxGradient()
{
    return mMaxGradient;
}

size_t KamadaKawaiReducer::FindMaxGradient::GetMaxGradientIdx()
{
    return mMaxGradientIdx;
}

KamadaKawaiReducer::IterateNewtonRaphson::IterateNewtonRaphson(
    double separationFactor, size_t molIdx,
    MolPtrVector &mols, double **dist, double **stiff
    ) :
    mFactor(separationFactor),
    mMolIdx(molIdx),
    mMols(mols),
    mDist(dist),
    mStiff(stiff),
    dEdxdx(0), dEdxdy(0),
    dEdydx(0), dEdydy(0),
    dEdx(0), dEdy(0)
{
    assert(mDist);
    assert(mStiff);
}

KamadaKawaiReducer::IterateNewtonRaphson::IterateNewtonRaphson(
    IterateNewtonRaphson &toSplit, tbb::split
    ) :
    mFactor(toSplit.mFactor),
    mMolIdx(toSplit.mMolIdx),
    mMols(toSplit.mMols),
    mDist(toSplit.mDist),
    mStiff(toSplit.mStiff),
    dEdxdx(0), dEdxdy(0),
    dEdydx(0), dEdydy(0),
    dEdx(0), dEdy(0)
{
}

void KamadaKawaiReducer::IterateNewtonRaphson::operator()(
    const tbb::blocked_range<size_t> &r)
{
    // Calculate Jacobi matrix (first-order partial derivatives).
    for (size_t i = r.begin(); i != r.end(); ++i) {
        double dx = mMols[mMolIdx]->posX - mMols[i]->posX;
        double dy = mMols[mMolIdx]->posY - mMols[i]->posY;
        double euclidDistance = sqrt(dx * dx + dy * dy);
        double cubedDistance = euclidDistance * euclidDistance * euclidDistance;

        if (cubedDistance > 0) {
            dEdxdx += mStiff[mMolIdx][i] *
                (1 - (mFactor * mDist[mMolIdx][i] * dy * dy) / cubedDistance);
            dEdxdy += mStiff[mMolIdx][i] *
                (mFactor * mDist[mMolIdx][i] * dx * dy) / cubedDistance;
            dEdydx = dEdxdy;
            dEdydy += mStiff[mMolIdx][i] *
                (1 - (mFactor * mDist[mMolIdx][i] * dx * dx) / cubedDistance);
            dEdx += mStiff[mMolIdx][i] *
                (dx - (mFactor * mDist[mMolIdx][i] * dx) / euclidDistance);
            dEdy += mStiff[mMolIdx][i] *
                (dy - (mFactor * mDist[mMolIdx][i] * dy) / euclidDistance);
        }
    }
}

void KamadaKawaiReducer::IterateNewtonRaphson::join(IterateNewtonRaphson &toJoin)
{
    dEdxdx += toJoin.dEdxdx;
    dEdxdy += toJoin.dEdxdy;
    dEdydx += toJoin.dEdydx;
    dEdydy += toJoin.dEdydy;
    dEdx += toJoin.dEdx;
    dEdy += toJoin.dEdy;
}

double KamadaKawaiReducer::IterateNewtonRaphson::GetGradient()
{
    return sqrt(dEdx * dEdx + dEdy * dEdy);
}

void KamadaKawaiReducer::IterateNewtonRaphson::UpdateCoordinates(
    double &x, double &y)
{
    double _dEdx = -dEdx;
    double _dEdy = -dEdy;

    // Cramer's rule
    double denom = dEdxdx * dEdydy - dEdxdy * dEdydx;
    if (denom != 0) {
        x += (_dEdx * dEdydy - dEdxdy * _dEdy) / denom;
        y += (dEdxdx * _dEdy - _dEdx * dEdydx) / denom;
    }
}

void KamadaKawaiReducer::Reduce(
        MolPtrVector &mols,
        FingerprintSelector fingerprintSelector,
        SimCoeffSelector simCoeffSelector,
        tbb::task_group_context &tbbCtx)
{
    /*
     For theoretical background and explanation of the algorithm, see
     http://www.koupy.net/download/GraphRec-thesis.pdf, section 2.1.2
     */

    MeasureStage measureStage;
    SimCoefCalculator calc(simCoeffSelector, fingerprintSelector);

    RandomizeCoordinates randomizeCoordinates;
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<MolPtrVector::iterator>(mols.begin(), mols.end()),
            randomizeCoordinates, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("RandomizeCoordinates");
    }

    std::vector<Fingerprint *> fingerprints;
    fingerprints.resize(mols.size(), NULL);
    CalculateFingerprints calculateFingerprints(calc, mols, fingerprints);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mols.size()),
            calculateFingerprints, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateFingerprints");
    }

    double **dist;
    dist = new double *[mols.size()];
    for (size_t i = 0; i < mols.size(); ++i) {
        dist[i] = new double[mols.size()];
    }
    CalculateDistances calculateDistances(calc, mols, fingerprints, dist);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range2d<size_t>(0, mols.size(), 0, mols.size()),
            calculateDistances, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateDistances");
    }
    for (size_t i = 0; i < mols.size(); ++i) {
        delete fingerprints[i];
    }
    fingerprints.clear();

    double **temp;
    temp = new double *[mols.size()];
    for (size_t i = 0; i < mols.size(); ++i) {
        temp[i] = new double[mols.size()];
    }
    // Find all-pairs shortest paths.
    for (size_t k = 0; k < mols.size(); ++k) {

        if (!Cancelled(tbbCtx)) {
            IterateFloydWarshall iterateFloydWarshall(k, dist, temp);
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, mols.size(), 0, mols.size()),
                iterateFloydWarshall, tbb::auto_partitioner(), tbbCtx);
            double **swap = temp;
            temp = dist;
            dist = swap;
        }

        if (Cancelled(tbbCtx)) {
            break;
        }
    }
    for (size_t i = 0; i < mols.size(); ++i) {
        delete[] temp[i];
    }
    delete[] temp;
    if (!Cancelled(tbbCtx)) {
        measureStage.ReportAndReset("FloydWarshall");
    }

    double **stiff;
    stiff = new double *[mols.size()];
    for (size_t i = 0; i < mols.size(); ++i) {
        stiff[i] = new double[mols.size()];
    }
    CalculateStiffness calculateStiffness(KAMADAKAWAIREDUCER_STIFFNESS, dist, stiff);
    if (!Cancelled(tbbCtx)) {
        tbb::parallel_for(
            tbb::blocked_range2d<size_t>(0, mols.size(), 0, mols.size()),
            calculateStiffness, tbb::auto_partitioner(), tbbCtx);
        measureStage.ReportAndReset("CalculateStiffness");
    }

    double currentGradient = DBL_MAX;
    size_t currentMolIdx = 0;
    size_t outCycles = 0;
    do {

        if (currentGradient != DBL_MAX) {
            // Find zero of nonlinear system of equations.
            size_t inCycles = 0;
            while (currentGradient > KAMADAKAWAIREDUCER_LOCALEPSILON) {

                IterateNewtonRaphson iterateNewtonRaphson(
                    KAMADAKAWAIREDUCER_SEPARATION, currentMolIdx, mols, dist, stiff);
                if (!Cancelled(tbbCtx)) {
                    tbb::parallel_reduce(
                        tbb::blocked_range<size_t>(0, mols.size()),
                        iterateNewtonRaphson, tbb::auto_partitioner(), tbbCtx);
                    //measureStage.ReportAndReset("IterateNewtonRaphson");
                }
                currentGradient = iterateNewtonRaphson.GetGradient();
                iterateNewtonRaphson.UpdateCoordinates(
                    mols[currentMolIdx]->posX, mols[currentMolIdx]->posY);

                if (++inCycles > KAMADAKAWAIREDUCER_MAXLOCALITER) {
                    mols[currentMolIdx]->posX = SynchRand::GetRandomNumber(0, 800);
                    mols[currentMolIdx]->posY = SynchRand::GetRandomNumber(0, 800);
                    break;
                }

                if (Cancelled(tbbCtx)) {
                    break;
                }
            }
        }
        //if (!Cancelled(tbbCtx)) {
        //    measureStage.ReportAndReset("IterateNewtonRaphson");
        //}

        // Determine molecule for next round.
        FindMaxGradient findMaxGradient(
            KAMADAKAWAIREDUCER_SEPARATION, mols, dist, stiff);
        if (!Cancelled(tbbCtx)) {
            tbb::parallel_reduce(
                tbb::blocked_range<size_t>(0, mols.size()),
                findMaxGradient, tbb::auto_partitioner(), tbbCtx);
            //measureStage.ReportAndReset("FindMaxGradient");
        }
        currentGradient = findMaxGradient.GetMaxGradient();
        currentMolIdx = findMaxGradient.GetMaxGradientIdx();

        if (Cancelled(tbbCtx) ||
                (++outCycles > mols.size() * KAMADAKAWAIREDUCER_GLOBALITERRATIO)) {
            break;
        }
    } while (currentGradient > KAMADAKAWAIREDUCER_GLOBALEPSILON);
    if (!Cancelled(tbbCtx)) {
        measureStage.ReportAndReset("CalculateCoordinates");
    }

    for (size_t i = 0; i < mols.size(); ++i) {
        delete[] stiff[i];
    }
    delete[] stiff;
    for (size_t i = 0; i < mols.size(); ++i) {
        delete[] dist[i];
    }
    delete[] dist;
}

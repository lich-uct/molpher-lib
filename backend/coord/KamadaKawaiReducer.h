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

#pragma once

#include <string>
#include <vector>
#include <ctime>

#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>

#include <tbb/task.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_reduce.h>

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "DimensionReducer.h"
#include "chem/simCoefStrategy/SimCoefStrategy.h"
#include "chem/SimCoefCalculator.hpp"

#ifndef KAMADAKAWAIREDUCER_REPORTING
#define KAMADAKAWAIREDUCER_REPORTING 1
#endif

#ifndef KAMADAKAWAIREDUCER_STIFFNESS
#define KAMADAKAWAIREDUCER_STIFFNESS 40.0
#endif

#ifndef KAMADAKAWAIREDUCER_SEPARATION
#define KAMADAKAWAIREDUCER_SEPARATION 4.0
#endif

#ifndef KAMADAKAWAIREDUCER_GLOBALEPSILON
#define KAMADAKAWAIREDUCER_GLOBALEPSILON 5.0
#endif

#ifndef KAMADAKAWAIREDUCER_LOCALEPSILON
#define KAMADAKAWAIREDUCER_LOCALEPSILON 0.5
#endif

#ifndef KAMADAKAWAIREDUCER_GLOBALITERRATIO
#define KAMADAKAWAIREDUCER_GLOBALITERRATIO 100
#endif

#ifndef KAMADAKAWAIREDUCER_MAXLOCALITER
#define KAMADAKAWAIREDUCER_MAXLOCALITER 100
#endif

class KamadaKawaiReducer : public DimensionReducer
{
public:
    KamadaKawaiReducer();
    virtual ~KamadaKawaiReducer();

    virtual void Reduce(
        MolPtrVector &mols,
        FingerprintSelector fingerprintSelector,
        SimCoeffSelector simCoeffSelector,
        tbb::task_group_context &tbbCtx);

protected:
    class RandomizeCoordinates
    {
    public:
        void operator()(
            const tbb::blocked_range<MolPtrVector::iterator> &mols) const;
    };

    class CalculateFingerprints
    {
    public:
        CalculateFingerprints(SimCoefCalculator &calc,
            MolPtrVector &mols, std::vector<Fingerprint *> &fingerprints);
        void operator()(const tbb::blocked_range<size_t> &r) const;

    private:
        SimCoefCalculator &mCalc;
        MolPtrVector &mMols;
        std::vector<Fingerprint *> &mFingerprints;
    };

    class CalculateDistances
    {
    public:
        CalculateDistances(SimCoefCalculator &calc, MolPtrVector &mols,
            std::vector<Fingerprint *> &fingerprints, double **dist);
        void operator()(const tbb::blocked_range2d<size_t> &r) const;

    private:
        SimCoefCalculator &mCalc;
        MolPtrVector &mMols;
        std::vector<Fingerprint *> &mFingerprints;
        double **mDist;
    };

    class IterateFloydWarshall
    {
    public:
        IterateFloydWarshall(size_t iter, double **last, double **current);
        void operator()(const tbb::blocked_range2d<size_t> &r) const;

    private:
        size_t mIter;
        double **mLast;
        double **mCurrent;
    };

    class CalculateStiffness
    {
    public:
        CalculateStiffness(double stiffnessFactor, double **dist, double **stiff);
        void operator()(const tbb::blocked_range2d<size_t> &r) const;

    private:
        double mFactor;
        double **mDist;
        double **mStiff;
    };

    class FindMaxGradient
    {
    public:
        FindMaxGradient(double separationFactor,
            MolPtrVector &mols, double **dist, double **stiff);
        FindMaxGradient(FindMaxGradient &toSplit, tbb::split);
        void operator()(const tbb::blocked_range<size_t> &r);
        void join(FindMaxGradient &toJoin);
        double GetMaxGradient();
        size_t GetMaxGradientIdx();

    private:
        double mFactor;
        double mMaxGradient;
        size_t mMaxGradientIdx;
        MolPtrVector &mMols;
        double **mDist;
        double **mStiff;
    };

    class IterateNewtonRaphson
    {
    public:
        IterateNewtonRaphson(double separationFactor, size_t molIdx,
            MolPtrVector &mols, double **dist, double **stiff);
        IterateNewtonRaphson(IterateNewtonRaphson &toSplit, tbb::split);
        void operator()(const tbb::blocked_range<size_t> &r);
        void join(IterateNewtonRaphson &toJoin);
        double GetGradient();
        void UpdateCoordinates(double &x, double &y);

    private:
        double mFactor;
        size_t mMolIdx;
        MolPtrVector &mMols;
        double **mDist;
        double **mStiff;

        double dEdxdx; // dE/dxdx
        double dEdxdy; // dE/dxdy
        double dEdydx; // dE/dydx
        double dEdydy; // dE/dydy
        double dEdx; // dE/dx
        double dEdy; // dE/dy
    };

    class MeasureStage
    {
    public:
        MeasureStage();
        void ReportAndReset(const std::string &stage);

    private:
        clock_t mTimestamp;
    };

    bool Cancelled(tbb::task_group_context &ctx);
};

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
#include <string>
#include <sstream>
#include <ctime>
#include <vector>

#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_exception.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"
#include "coord/ReducerFactory.h"
#include "chem/morphing/Morphing.hpp"
#include "NeighborhoodTaskQueue.h"
#include "NeighborhoodGenerator.h"

NeighborhoodGenerator::NeighborhoodGenerator(
    tbb::task_group_context *tbbCtx, NeighborhoodTaskQueue *queue, int threadCnt
    ) :
    mTbbCtx(tbbCtx),
    mQueue(queue),
    mThreadCnt(threadCnt)
{
}

NeighborhoodGenerator::~NeighborhoodGenerator()
{
}

bool NeighborhoodGenerator::Cancelled()
{
    return mTbbCtx->is_group_execution_cancelled();
}

void NeighborhoodCollector(MolpherMolecule *morph, void *functor)
{
    NeighborhoodGenerator::CollectMorph *collect =
        (NeighborhoodGenerator::CollectMorph *) functor;
    (*collect)(*morph);
}

NeighborhoodGenerator::CollectMorph::CollectMorph() :
    mMorphSet(false)
{
}

void NeighborhoodGenerator::CollectMorph::operator()(const MolpherMolecule &morph)
{
    assert(mMorphSet == false);
    mMorph = morph;
    mMorphSet = true;
}

bool NeighborhoodGenerator::CollectMorph::WithdrawMorph(MolpherMolecule &morph)
{
    if (mMorphSet) {
        mMorphSet = false;
        morph = mMorph;
        return true;
    } else {
        return false;
    }
}

NeighborhoodGenerator::GenerateNeighborhood::GenerateNeighborhood(
    NeighborhoodTask &task, MoleculeVector &neighborhood,
    tbb::task_group_context *tbbCtx
    ) :
    mTask(task),
    mNeighborhood(neighborhood),
    mTbbCtx(tbbCtx)
{
}

void NeighborhoodGenerator::GenerateNeighborhood::operator()(
    const tbb::blocked_range<size_t> &r) const
{
    CollectMorph collectMorph;
    MolpherMolecule generator;
    MolpherMolecule neighbor;

    std::vector<ChemOperSelector> chemOperSelectors;
    chemOperSelectors.resize(mTask.chemOperSelectors.size(), (ChemOperSelector) 0);
    for (size_t i = 0; i < mTask.chemOperSelectors.size(); ++i) {
        chemOperSelectors[i] = (ChemOperSelector) mTask.chemOperSelectors[i];
    }

    std::vector<MolpherMolecule> emptyDecoys;

    for (size_t attempt = r.begin(); attempt != r.end(); ++attempt) {

        int depth = SynchRand::GetRandomNumber(1, mTask.maxDepth);
        generator = mTask.origin;

        for (int d = 0; d < depth; ++d) {
            while (!collectMorph.WithdrawMorph(neighbor)) {
                GenerateMorphs(
                    generator,
                    1,
                    (FingerprintSelector) mTask.fingerprintSelector,
                    (SimCoeffSelector) mTask.simCoeffSelector,
                    chemOperSelectors,
                    mTask.origin,
                    emptyDecoys,
                    *mTbbCtx,
                    &collectMorph,
                    NeighborhoodCollector);
            }
            generator = neighbor;
        }

        bool withinNeighborhood = (neighbor.distToTarget <= mTask.maxDistance);
        if (withinNeighborhood) {
            mNeighborhood.push_back(neighbor);
        }
    }
}

void NeighborhoodGenerator::operator()()
{
    SynchCout(std::string("NeighborhoodGenerator thread started."));

    tbb::task_scheduler_init scheduler;
    if (mThreadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(mThreadCnt);
    }

    while (true) {

        NeighborhoodTask task;
        NeighborhoodTaskResult result;

        if (!mQueue->Pop(task)) {
            break; // Thread termination.
        }

        try {

            /* TODO MPI
             broadcast task
             scatter attempts over cluster
            */

            MoleculeVector neighborhood;
            GenerateNeighborhood generateNeighborhood(task, neighborhood, mTbbCtx);
            if (!Cancelled() && !task.origin.smile.empty()) {

                clock_t start = std::clock();
                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, task.attemptCount),
                    generateNeighborhood, tbb::auto_partitioner(), *mTbbCtx);
                clock_t finish = std::clock();

#if NEIGHBORHOODGENERATOR_REPORTING == 1
                std::ostringstream stream;
                stream << boost::posix_time::to_iso_string(task.taskTimestamp) <<
                    ": " << "GenerateNeighborhood consumed " <<
                    finish - start << " msec.";
                SynchCout(stream.str());
#endif
            }

            /* TODO MPI
             convert neighborhood to std::vector
             gather neighborhood over cluster
            */

            if (!Cancelled()) {
                result.taskTimestamp = task.taskTimestamp;
                result.origin = task.origin;
                result.reducedNeighborhood.insert(
                    result.reducedNeighborhood.end(),
                    neighborhood.begin(), neighborhood.end());
                result.reducedContext = task.context;
            }

            if (!Cancelled()) {
                clock_t start = std::clock();

                DimensionReducer::MolPtrVector molsToReduce;
                molsToReduce.reserve(result.reducedNeighborhood.size() +
                    result.reducedContext.size() + 1);
                std::vector<MolpherMolecule>::iterator itNeighborhood;
                for (itNeighborhood = result.reducedNeighborhood.begin();
                        itNeighborhood != result.reducedNeighborhood.end();
                        itNeighborhood++) {
                    molsToReduce.push_back(&(*itNeighborhood));
                }
                std::vector<MolpherMolecule>::iterator itContext;
                for (itContext = result.reducedContext.begin();
                        itContext != result.reducedContext.end(); itContext++) {
                    molsToReduce.push_back(&(*itContext));
                }
                if (!result.origin.smile.empty()) {
                    molsToReduce.push_back(&result.origin);
                }

                DimensionReducer *reducer =
                    ReducerFactory::Create((DimRedSelector) task.dimRedSelector);
                reducer->Reduce(molsToReduce,
                    (FingerprintSelector) task.fingerprintSelector,
                    (SimCoeffSelector) task.simCoeffSelector, *mTbbCtx);
                ReducerFactory::Recycle(reducer);

                clock_t finish = std::clock();

#if NEIGHBORHOODGENERATOR_REPORTING == 1
                std::ostringstream stream;
                stream << boost::posix_time::to_iso_string(task.taskTimestamp) <<
                    ": " << "DimensionReduction consumed " <<
                    finish - start << " msec.";
                SynchCout(stream.str());
#endif
            }

        } catch (tbb::tbb_exception &exc) {
            SynchCout(std::string(exc.what()));
            mTbbCtx->cancel_group_execution(); // Invalidate result.
        }

        mQueue->CommitTaskResult(result);

    }

    SynchCout(std::string("NeighborhoodGenerator thread terminated."));
}

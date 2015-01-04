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

#include <tbb/task.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

#include "NeighborhoodTask.h"

#ifndef NEIGHBORHOODGENERATOR_REPORTING
#define NEIGHBORHOODGENERATOR_REPORTING 1
#endif

class NeighborhoodTaskQueue;

class NeighborhoodGenerator
{
public:
    NeighborhoodGenerator(tbb::task_group_context *tbbCtx,
        NeighborhoodTaskQueue *queue, int threadCnt = 0);
    ~NeighborhoodGenerator();

    void operator()();

protected:
    typedef tbb::concurrent_vector<MolpherMolecule> MoleculeVector;

    friend void NeighborhoodCollector(MolpherMolecule *morph, void *functor);

    class CollectMorph
    {
    public:
        CollectMorph();
        void operator()(const MolpherMolecule &morph);
        bool WithdrawMorph(MolpherMolecule &morph);

    private:
        bool mMorphSet;
        MolpherMolecule mMorph;
    };

    class GenerateNeighborhood
    {
    public:
        GenerateNeighborhood(NeighborhoodTask &task,
            MoleculeVector &neighborhood, tbb::task_group_context *tbbCtx);
        void operator()(const tbb::blocked_range<size_t> &r) const;

    private:
        NeighborhoodTask &mTask;
        MoleculeVector &mNeighborhood;
        tbb::task_group_context *mTbbCtx;
    };

    bool Cancelled();

private:
    tbb::task_group_context *mTbbCtx;
    NeighborhoodTaskQueue *mQueue;
    int mThreadCnt;
};

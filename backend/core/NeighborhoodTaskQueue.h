/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Vladimir Fiklik

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

#include <deque>

#include <tbb/task.h>

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>

#include "NeighborhoodTask.h"

class BackendCommunicator;

class NeighborhoodTaskQueue
{
public:
    // Functions called by backend main thread.
    NeighborhoodTaskQueue(tbb::task_group_context *neighborhoodGeneratorStopper);
    ~NeighborhoodTaskQueue();
    void SetCommunicator(BackendCommunicator *comm);
    void Halt();

    // Functions called by neighborhood generator top-level thread.
    bool Pop(NeighborhoodTask &task);
    void CommitTaskResult(NeighborhoodTaskResult &res);

    // Functions called by communicator thread.
    void Push(NeighborhoodTask &task);
    void SkipNeighborhoodTask(boost::posix_time::ptime timestamp);

protected:
    // Functions called only internally. Assumes proper synchronization by caller.
    void PublishTaskResult(NeighborhoodTaskResult &res);

private:
    typedef boost::mutex Guard;
    typedef boost::unique_lock<Guard> Lock;

    BackendCommunicator *mCommunicator; // Provides publishing functionality.
    tbb::task_group_context *mNeighborhoodGeneratorStopper; // Flushes current task.
    bool mHalted; // Used for neighborhood generator thread termination.
    boost::condition_variable mTaskReadyCondition; // Wakes sleeping neighborhood generator.
    Guard mNeighborhoodTaskQueueGuard; // Serializes thread access to most of the methods.

    boost::posix_time::ptime mCurrentTaskTimestamp;
    std::deque<NeighborhoodTask> mTaskQueue;
};

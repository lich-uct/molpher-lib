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

#include "inout.h"
#include "BackendCommunicator.h"
#include "NeighborhoodTaskQueue.h"

NeighborhoodTaskQueue::NeighborhoodTaskQueue(
    tbb::task_group_context *neighborhoodGeneratorStopper
    ) :
    mHalted(false),
    mNeighborhoodGeneratorStopper(neighborhoodGeneratorStopper),
    mCommunicator(0),
    mCurrentTaskTimestamp(boost::date_time::min_date_time)
{
    SynchCout(std::string("NeighborhoodTaskQueue initialized."));
}

void NeighborhoodTaskQueue::SetCommunicator(BackendCommunicator *comm)
{
    mCommunicator = comm;
}

void NeighborhoodTaskQueue::Halt()
{
    Lock lock(mNeighborhoodTaskQueueGuard);
    mHalted = true;
    mTaskQueue.clear(); // Flush the queue, so there is nothing else to pop.
    mNeighborhoodGeneratorStopper->cancel_group_execution();
    mTaskReadyCondition.notify_all(); // In case the neighborhood generator is sleeping.

    SynchCout(std::string("NeighborhoodTaskQueue halted."));
}

NeighborhoodTaskQueue::~NeighborhoodTaskQueue()
{
}

bool NeighborhoodTaskQueue::Pop(NeighborhoodTask &task)
{
    Lock lock(mNeighborhoodTaskQueueGuard);
    while(mTaskQueue.empty()) {
        if (mHalted) {
            return false; // Neighborhood generator thread will terminate.
        } else {
            mTaskReadyCondition.wait(lock); // Yields lock until signalled.
        }
    }
    task = mTaskQueue.front(); // Copy the task.
    mCurrentTaskTimestamp = task.taskTimestamp; // Remember the current task's timestamp
    mTaskQueue.pop_front();
    return true;
}

void NeighborhoodTaskQueue::CommitTaskResult(NeighborhoodTaskResult &res)
{
    Lock lock(mNeighborhoodTaskQueueGuard);

    mCurrentTaskTimestamp = boost::date_time::min_date_time;
    if (mNeighborhoodGeneratorStopper->is_group_execution_cancelled()) {
        mNeighborhoodGeneratorStopper->reset();
    } else {
        PublishTaskResult(res);
    }
}

void NeighborhoodTaskQueue::Push(NeighborhoodTask &task)
{
    Lock lock(mNeighborhoodTaskQueueGuard);
    if (task.IsValid()) { // Prevents backend crash (should be ensured by frontend).
        mTaskQueue.push_back(task);
    }
    lock.unlock(); // Unlock to prevent deadlock when signalling the condition.
    mTaskReadyCondition.notify_all();
}

void NeighborhoodTaskQueue::SkipNeighborhoodTask(boost::posix_time::ptime timestamp)
{
    Lock lock(mNeighborhoodTaskQueueGuard);
    if (timestamp == mCurrentTaskTimestamp) {
        // If the timestamp matches the current task, stop the current task
        mNeighborhoodGeneratorStopper->cancel_group_execution();
    } else {
        // Else search through the queue and remove the appropriate task
        std::deque<NeighborhoodTask>::iterator it = mTaskQueue.begin();
        while (it != mTaskQueue.end()) {
            if ((*it).taskTimestamp == timestamp) {
                mTaskQueue.erase(it);
                break;
            }
        }
    }
}

void NeighborhoodTaskQueue::PublishTaskResult(NeighborhoodTaskResult &res)
{
    // Already locked by caller.
    if (mCommunicator) {
        mCommunicator->PublishNeighborhoodTaskResult(res);
    }
}

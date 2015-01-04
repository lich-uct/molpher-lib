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

#include <cassert>
#include <map>
#include <deque>

#include <boost/cstdint.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include "IterationSnapshot.h"

struct JobGroup
{
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & mJobMap & mLiveJobQueue & mSleepingJobQueue & mFinishedJobQueue;
    }

    typedef std::map<boost::uint32_t, IterationSnapshot> JobMap;
    typedef std::deque<boost::uint32_t> JobQueue;

    /*
     job states:
     - live jobs
         - are one by one consumed by path finder
         - first job in live queue is implicitly considered as running
     - sleeping jobs
         - not queued, but can be woken up to become live job
         - job can be put to sleep by request from frontend or by path finder
           if such job consumed all iterations or seconds specified inside its
           last executed iteration without finding the target (those params
           could be however later edited from frontend and job could be
           consequently revived)
     - terminated jobs
         - successfully finished jobs (target found, path established)
     */

    JobMap mJobMap; // maps job ID to the most recent iteration of the job
    JobQueue mLiveJobQueue; // stores order in which jobs are being computed
    JobQueue mSleepingJobQueue; // stores order in which jobs were deactivated
    JobQueue mFinishedJobQueue; // stores order in which jobs were finished

    void GetJob(unsigned int idx, IterationSnapshot &snp)
    {
        if (idx < mLiveJobQueue.size()) {
            snp = mJobMap[mLiveJobQueue[idx]];
        } else if (idx < mLiveJobQueue.size() + mSleepingJobQueue.size()) {
            snp = mJobMap[mSleepingJobQueue[idx - mLiveJobQueue.size()]];
        } else if (idx < mLiveJobQueue.size() + mSleepingJobQueue.size() + mFinishedJobQueue.size()) {
            snp = mJobMap[mFinishedJobQueue[idx - mLiveJobQueue.size() - mSleepingJobQueue.size()]];
        } else {
            assert(false);
        }
    }

    unsigned int GetIndex(boost::uint32_t jobId)
    {
        for (size_t i = 0; i < mLiveJobQueue.size(); ++i) {
            if (mLiveJobQueue[i] == jobId) {
                return i;
            }
        }

        for (size_t i = 0; i < mSleepingJobQueue.size(); ++i) {
            if (mSleepingJobQueue[i] == jobId) {
                return mLiveJobQueue.size() + i;
            }
        }

        for (size_t i = 0; i < mFinishedJobQueue.size(); ++i) {
            if (mFinishedJobQueue[i] == jobId) {
                return mLiveJobQueue.size() + mSleepingJobQueue.size() + i;
            }
        }

        assert(false);
        return -1;
    }

    bool IsJobInQueue(JobGroup::JobQueue &queue, boost::uint32_t jobId)
    {
        for (size_t i = 0; i < queue.size(); ++i) {
            if (queue[i] == jobId) {
                return true;
            }
        }
        return false;
    }

    bool IsLive(boost::uint32_t jobId)
    {
        return IsJobInQueue(mLiveJobQueue, jobId);
    }

    bool IsSleeping(boost::uint32_t jobId)
    {
        return IsJobInQueue(mSleepingJobQueue, jobId);
    }

    bool IsFinished(boost::uint32_t jobId)
    {
        return IsJobInQueue(mFinishedJobQueue, jobId);
    }
};

BOOST_CLASS_IMPLEMENTATION(JobGroup, object_serializable) // turn off versioning
BOOST_CLASS_TRACKING(JobGroup, track_never) // turn off tracking

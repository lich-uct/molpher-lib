/*
 Copyright (c) 2012 Vladimir Fiklik
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

#include <boost/cstdint.hpp>

#include <RCF/RCF.hpp>
#include <RCF/PublishingService.hpp>
#include <RCF/FilterService.hpp>
#include <RCF/ZlibCompressionFilter.hpp>

#include "MolpherParam.h"
#include "MolpherMolecule.h"
#include "IterationSnapshot.h"
#include "JobGroup.h"
#include "NeighborhoodTask.h"
#include "core/NeighborhoodTaskQueue.h"
#include "core/JobManager.h"

class JobManager;
class BackendCommunicator;

class FrontendConnectedHandler
{
public:
    FrontendConnectedHandler(BackendCommunicator *communicator);
    void operator()(RCF::RcfSession &session, const std::string &iface);

private:
    BackendCommunicator *mCommunicator;
};

class FrontendDisconnectedHandler
{
public:
    FrontendDisconnectedHandler(BackendCommunicator *communicator);
    void operator()(RCF::RcfSession &session, const std::string &iface);

private:
    BackendCommunicator *mCommunicator;
};

class BackendCommunicator
{
public:
    BackendCommunicator(JobManager *jobManager, NeighborhoodTaskQueue *taskQueue);
    ~BackendCommunicator();
    void Halt();

    // Following methods must exactly match the BackendIfc.
    void InitClient(std::string &backendId);
    boost::uint32_t CreateJob(IterationSnapshot snp, std::string password);
    void WakeJob(boost::uint32_t jobId, std::string password);
    void SleepJob(boost::uint32_t jobId, std::string password);
    void RemoveJob(boost::uint32_t jobId, std::string password);
    void ChangeJobOrder(boost::uint32_t jobId, boost::int32_t queuePosDiff, std::string password);
    bool ValidateJobPassword(boost::uint32_t jobId, std::string password);
    IterationSnapshot GetJobHistory(boost::uint32_t jobId, boost::uint32_t iterIdx, bool &loaded);
    bool SetFingerprintSelector(boost::uint32_t jobId, boost::int32_t selector, std::string password);
    bool SetSimCoeffSelector(boost::uint32_t jobId, boost::int32_t selector, std::string password);
    bool SetDimRedSelector(boost::uint32_t jobId, boost::int32_t selector, std::string password);
    bool SetChemOperSelectors(boost::uint32_t jobId, std::vector<boost::int32_t> selectors,
        std::string password);
    bool SetParams(boost::uint32_t jobId, MolpherParam params, std::string password);
    bool SetDecoys(boost::uint32_t jobId, std::vector<MolpherMolecule> decoys, std::string password);
    bool AddPruned(boost::uint32_t jobId, std::vector<MolpherMolecule> pruned, std::string password);
    void EnqueueNeighborhoodTask(NeighborhoodTask task);
    void SkipNeighborhoodTask(boost::posix_time::ptime timestamp);

    void PublishJobs(JobGroup &jobs);
    void PublishIteration(IterationSnapshot &snp);
    void PublishNeighborhoodTaskResult(NeighborhoodTaskResult &res);

private:
    RCF::RcfServer mPublisher;
    RCF::RcfServer mListener;
    RCF::PublishingServicePtr mPubSvc;
    RCF::FilterServicePtr mFltSvc;
    RCF::FilterPtr mCompressFlt;
    FrontendConnectedHandler mConnectHandler;
    FrontendDisconnectedHandler mDisconnectHandler;

    JobManager *mJobManager;
    NeighborhoodTaskQueue *mTaskQueue;
};

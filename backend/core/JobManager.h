/*
 Copyright (c) 2012 Petr Koupy
 Copyrighy (c) 2012 Vladimir Fiklik

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

#include <tbb/task.h>

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition_variable.hpp>

#include "MolpherParam.h"
#include "MolpherMolecule.h"
#include "PathFinderContext.h"
#include "IterationSnapshot.h"
#include "JobGroup.h"

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "dimred_selectors.h"
#include "chemoper_selectors.h"

class BackendCommunicator;

class JobManager
{
public:
    // Functions called by backend main thread.
    JobManager(tbb::task_group_context *pathFinderStopper,
        std::string &storagePath, std::string &jobListFile, bool interactive);
    ~JobManager();
    void AddJobFromFile(const std::string &jobFile);
    void SetCommunicator(BackendCommunicator *comm);
    void Halt();

    // Functions called by path finder top-level thread.
    bool GetJob(PathFinderContext &ctx);
    bool CommitIteration(PathFinderContext &ctx, bool canContinue, bool pathFound);
    bool GetFingerprintSelector(FingerprintSelector &selector);
    bool GetSimCoeffSelector(SimCoeffSelector &selector);
    bool GetDimRedSelector(DimRedSelector &selector);
    bool GetChemOperSelectors(std::vector<ChemOperSelector> &selectors);
    bool GetParams(MolpherParam &params);
    bool GetDecoys(std::vector<MolpherMolecule> &decoys);
    void GetPruned(std::vector<MolpherMolecule> &pruned);

    // Functions called by communicator thread.
    void OnConnect(std::string &backendId);
    JobId CreateJob(IterationSnapshot &snp, std::string &password);
    void WakeJob(JobId jobId, std::string &password);
    void SleepJob(JobId jobId, std::string &password);
    void RemoveJob(JobId jobId, std::string &password);
    void ChangeJobOrder(JobId jobId, int queuePosDiff, std::string &password);
    bool ValidateJobPassword(JobId jobId, std::string &password);
    bool GetJobHistory(JobId jobId, IterIdx iterIdx, IterationSnapshot &snp);
    bool SetFingerprintSelector(JobId jobId, FingerprintSelector selector, std::string &password);
    bool SetSimCoeffSelector(JobId jobId, SimCoeffSelector selector, std::string &password);
    bool SetDimRedSelector(JobId jobId, DimRedSelector selector, std::string &password);
    bool SetChemOperSelectors(JobId jobId, std::vector<ChemOperSelector> &selectors, std::string &password);
    bool SetParams(JobId jobId, MolpherParam &params, std::string &password);
    bool SetDecoys(JobId jobId, std::vector<MolpherMolecule> &decoys, std::string &password);
    bool AddPruned(JobId jobId, std::vector<MolpherMolecule> &pruned, std::string &password);

protected:
    // Functions called only internally. Assumes proper synchronization by caller.
    void PublishJobs();
    void PublishIteration(IterationSnapshot &snp);
    bool VerifyPassword(JobId jobId, std::string &password);
    bool GetFirstLiveJob(IterationSnapshot &snp);
    void DeleteJob(JobId jobId);
    bool DeleteFromQueue(JobGroup::JobQueue &queue, JobId jobId);
    void ClearDeferredActions();
    bool IsJobRunning(JobId jobId);

private:
    typedef boost::mutex Guard; // Could be changed to recursive mutex if needed.
    typedef boost::unique_lock<Guard> Lock;

    BackendCommunicator *mCommunicator; // Provides publishing functionality.
    tbb::task_group_context *mPathFinderStopper; // Flushes current iteration.
    boost::condition_variable mJobReadyCondition; // Wakes sleeping path finder.
    bool mHalted; // Used for path finder thread termination.
    bool mInteractive; // Non-interactive mode terminates on empty queue.
    Guard mJobManagerGuard; // Serializes thread access to most of the methods.
    Guard mDeferredActionsGuard; // Protects just deferred action data.

    JobId mJobIdCounter; // Unique IDs during single execution of backend.
    JobGroup mJobs; // Actual job queues and descriptions.

    FingerprintSelector mDeferredFingerprintSelector;
    SimCoeffSelector mDeferredSimCoeffSelector;
    DimRedSelector mDeferredDimRedSelector;
    std::vector<ChemOperSelector> mDeferredChemOperSelectors;
    MolpherParam mDeferredMolpherParams;
    std::vector<MolpherMolecule> mDeferredDecoys;
    std::vector<MolpherMolecule> mDeferredPrunedAccumulator;

    bool mDeferredFingerprintSelectorIsSet;
    bool mDeferredSimCoeffSelectorIsSet;
    bool mDeferredDimRedSelectorIsSet;
    bool mDeferredChemOperSelectorsIsSet;
    bool mDeferredMolpherParamsIsSet;
    bool mDeferredDecoysIsSet;

    typedef std::map<JobId, std::string> PasswordMap;
    PasswordMap mPasswordMap;

    std::string mStorageDir; // Dir to store the results
    std::string mBackendId; // Identification timestamp
};

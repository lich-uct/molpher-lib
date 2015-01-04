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

#include "auxiliary/QtMocHack.h"

#include <string>
#include <vector>

#include <boost/cstdint.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include <RCF/RCF.hpp>
#include <RCF/TcpEndpoint.hpp>
#include <RCF/SubscriptionService.hpp>
#include <RCF/FilterService.hpp>
#include <RCF/ZlibCompressionFilter.hpp>
#include <RCF/ClientProgress.hpp>

#include <QtCore/QObject>
#include <QtCore/QString>

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "dimred_selectors.h"
#include "chemoper_selectors.h"

#include "MolpherParam.h"
#include "MolpherMolecule.h"
#include "meta/IterationSnapshotMeta.h"
#include "meta/JobGroupMeta.h"
#include "meta/NeighborhoodTaskMeta.h"

extern std::string gFrontendId;

class FrontendCommunicator;
extern FrontendCommunicator gCommunicator;

class BackendDisconnectedHandler
{
public:
    BackendDisconnectedHandler(FrontendCommunicator *communicator);
    void operator()(RCF::RcfSession &session);

private:
    FrontendCommunicator *mCommunicator;
};

class FrontendCommunicator : public QObject
{
    Q_OBJECT
public:
    bool Init(); // Must be called before calling any other method.
    void DeInit(); // Must be called before calling destructor.

    void BackendDisconnected();

    // Following methods must exactly match the FrontendIfc.
    void AcceptJobs(JobGroup jobs);
    void AcceptIteration(IterationSnapshot snp);
    void AcceptNeighborhoodTaskResult(NeighborhoodTaskResult res);

public slots:
    // Slots that should be connected to GUI by Qt::DirectConnection.
    void ConnectToBackend(const std::string &ip);
    void DisconnectFromBackend();
    void CreateJob(IterationSnapshot &snp, std::string &password, JobId &jobId);
    void WakeJob(JobId jobId, std::string &password);
    void SleepJob(JobId jobId, std::string &password);
    void RemoveJob(JobId jobId, std::string &password);
    void ChangeJobOrder(JobId jobId, int queuePosDiff, std::string &password);
    void ValidateJobPassword(JobId jobId, std::string &password, bool &isValid);
    void GetJobHistory(JobId jobId, IterIdx minIterIdx, IterIdx maxIterIdx, bool full, std::vector<IterationSnapshotProxy> &history);
    void SetFingerprintSelector(JobId jobId, FingerprintSelector selector, std::string &password);
    void SetSimCoeffSelector(JobId jobId, SimCoeffSelector selector, std::string &password);
    void SetDimRedSelector(JobId jobId, DimRedSelector selector, std::string &password);
    void SetChemOperSelectors(JobId jobId, std::vector<ChemOperSelector> &selectors,
        std::string &password);
    void SetParams(JobId jobId, MolpherParam &params, std::string &password);
    void SetDecoys(JobId jobId, std::vector<MolpherMolecule> &decoys, std::string &password);
    void AddPruned(JobId jobId, std::vector<MolpherMolecule> &pruned, std::string &password);
    void EnqueueNeighborhoodTask(NeighborhoodTask &task);
    void SkipNeighborhoodTask(boost::posix_time::ptime timestamp);

signals:
    // Signals that should be connected to GUI by Qt::QueuedConnection.
    void DisplayOnlineState(const QString &ip);
    void DisplayOfflineState();
    void UpdateJobs(const JobGroup &jobs);
    void VisualizeIteration(const IterationSnapshot &snp);
    void VisualizeIteration(const IterationSnapshotProxy &proxy);
    void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);
    // Signals that should be connected to GUI by Qt::DirectConnection.
    void SendHistory(const IterationSnapshotProxy &proxy);
    void ReportStatus(const QString &message, int timeout);

protected:
    static void OnProgress(RCF::ClientProgress::Action &action);

private:
    typedef boost::mutex Guard;
    typedef boost::unique_lock<Guard> Lock;

    RCF::RcfServer *mSubscriber;
    RCF::SubscriptionServicePtr mSubSvc;
    RCF::FilterServicePtr mFltSvc;
    RCF::FilterPtr mCompressFlt;
    BackendDisconnectedHandler *mDisconnectHandler;
    RCF::TcpEndpoint mBackendEndpoint;
    Guard mBackendIdGuard;
    std::string mBackendId;
    std::string mStorageDir;
};

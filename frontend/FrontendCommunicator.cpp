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

#include <cassert>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <QtGui/QApplication>

#include "inout.h"
#include "auxiliary/PasswordCache.h"
#include "molpher_interface.idl"
#include "FrontendCommunicator.h"

#define MESSAGE_SIZE (1024 * 1024 * 1024) // 1GB

#define REPORT_TIMEOUT 5000

std::string gFrontendId(
    boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time()));

FrontendCommunicator gCommunicator;

bool FrontendCommunicator::Init()
{
    mStorageDir = "Temp";
    mStorageDir += "/" + gFrontendId;
    try {
        boost::filesystem::create_directories(mStorageDir);
    } catch (boost::filesystem::filesystem_error &exc) {
        SynchCout(exc.what());
    }

    mSubscriber = new RCF::RcfServer(RCF::TcpEndpoint("localhost", -1));
    mSubSvc = RCF::SubscriptionServicePtr(new RCF::SubscriptionService());
    mFltSvc = RCF::FilterServicePtr(new RCF::FilterService());
    mCompressFlt = RCF::FilterPtr(new RCF::ZlibStatelessCompressionFilter());
    mDisconnectHandler = new BackendDisconnectedHandler(this);

    bool initialized = mSubscriber && mSubSvc.get() && mFltSvc.get() &&
        mCompressFlt.get() && mDisconnectHandler;

    if (initialized) {
        mFltSvc->addFilterFactory( RCF::FilterFactoryPtr(
            new RCF::ZlibStatelessCompressionFilterFactory()));

        mSubscriber->getServerTransport().setMaxMessageLength(MESSAGE_SIZE);
        mSubscriber->addService(mSubSvc);
        mSubscriber->addService(mFltSvc);
        mSubscriber->start();
    } else {
        this->DeInit();
    }

    return initialized;
}

void FrontendCommunicator::DeInit()
{
    try {
        boost::filesystem::remove_all(mStorageDir);
    } catch (boost::filesystem::filesystem_error &exc) {
        SynchCout(exc.what());
    }

    delete mDisconnectHandler;
    mSubSvc.reset();
    mFltSvc.reset();
    mCompressFlt.reset();
    delete mSubscriber;
}

void FrontendCommunicator::ConnectToBackend(const std::string &ip)
{
    mSubSvc->endSubscribe<FrontendIfc>(*this); // Disconnect if connected.

    Lock lock(mBackendIdGuard);
    RCF::SubscriptionPtr subscription;
    try {
        mBackendEndpoint = RCF::TcpEndpoint(ip, 50002);
        subscription = mSubSvc->beginSubscribe<FrontendIfc>(
            *this, RCF::TcpEndpoint(ip, 50001), *mDisconnectHandler);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }

    if (!subscription || !subscription->isConnected()) {
        SynchCout(std::string("Cannot connect to backend."));
    } else {
        SynchCout(std::string("Connected to backend."));
        emit DisplayOnlineState(QString::fromStdString(ip));
        try {
            RcfClient<BackendIfc>(mBackendEndpoint).InitClient(mBackendId);
            boost::filesystem::create_directories(mStorageDir + "/" + mBackendId);
        } catch (RCF::Exception &exc) {
            SynchCout(exc.getErrorString());
        } catch (boost::filesystem::filesystem_error &exc) {
            SynchCout(exc.what());
        }
    }
}

void FrontendCommunicator::DisconnectFromBackend()
{
    Lock lock(mBackendIdGuard);
    SynchCout(std::string("Disconnected from backend."));
    mSubSvc->endSubscribe<FrontendIfc>(*this);
    mBackendId.clear();
    emit DisplayOfflineState();
}

void FrontendCommunicator::BackendDisconnected()
{
    Lock lock(mBackendIdGuard);
    mBackendId.clear();
    SynchCout(std::string("Disconnected from backend."));
    emit DisplayOfflineState();
}

void FrontendCommunicator::AcceptJobs(JobGroup jobs)
{
    emit UpdateJobs(jobs);
}

void FrontendCommunicator::AcceptIteration(IterationSnapshot snp)
{
    Lock lock(mBackendIdGuard);
    if (mBackendId.empty()) {
        return;
    }
    std::string storage = mStorageDir + "/" + mBackendId;
    lock.unlock();
    try {
        boost::filesystem::create_directories(GenerateDirname(storage, snp.jobId));
    } catch (boost::filesystem::filesystem_error &exc) {
        SynchCout(exc.what());
    }
    WriteSnapshotToFile(
        GenerateFilename(storage, snp.jobId, snp.iterIdx), snp);
    emit VisualizeIteration(snp);
    IterationSnapshotProxy proxy(storage, snp.jobId, snp.iterIdx);
    emit VisualizeIteration(proxy);
}

void FrontendCommunicator::AcceptNeighborhoodTaskResult(NeighborhoodTaskResult res)
{
    emit VisualizeNeighborhoodTaskResult(res);
}

void FrontendCommunicator::CreateJob(
    IterationSnapshot &snp, std::string &password, JobId &jobId)
{
    try {
        RcfClient<BackendIfc> client(mBackendEndpoint);
        client.getClientStub().getTransport().setMaxMessageLength(MESSAGE_SIZE);
        client.getClientStub().setMessageFilters(mCompressFlt);
        jobId = client.CreateJob(snp, password);
        PasswordCache::CachePassword(jobId, password);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        jobId = -1;
    }
}

void FrontendCommunicator::WakeJob(JobId jobId, std::string &password)
{
    try {
        RcfClient<BackendIfc>(mBackendEndpoint).WakeJob(jobId, password);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }
}

void FrontendCommunicator::SleepJob(JobId jobId, std::string &password)
{
    try {
        RcfClient<BackendIfc>(mBackendEndpoint).SleepJob(jobId, password);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }
}

void FrontendCommunicator::RemoveJob(JobId jobId, std::string &password)
{
    try {
        RcfClient<BackendIfc>(mBackendEndpoint).RemoveJob(jobId, password);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }
}

void FrontendCommunicator::ChangeJobOrder(JobId jobId,
    int queuePosDiff, std::string &password)
{
    try {
        RcfClient<BackendIfc>(mBackendEndpoint).ChangeJobOrder(jobId, queuePosDiff, password);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }
}

void FrontendCommunicator::ValidateJobPassword(JobId jobId,
    std::string &password, bool &isValid)
{
    try {
        isValid = RcfClient<BackendIfc>(mBackendEndpoint).ValidateJobPassword(jobId, password);
        if (isValid) {
            PasswordCache::CachePassword(jobId, password);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        isValid = false;
    }
}

void FrontendCommunicator::OnProgress(RCF::ClientProgress::Action &action)
{
    QCoreApplication::processEvents(); // Keep GUI responsive.
    action = RCF::ClientProgress::Continue; // Continue the call.
}

void FrontendCommunicator::GetJobHistory(JobId jobId, IterIdx minIterIdx,
    IterIdx maxIterIdx, bool full, std::vector<IterationSnapshotProxy> &history)
{
    Lock lock(mBackendIdGuard);
    if (mBackendId.empty()) {
        return;
    }
    std::string storage = mStorageDir + "/" + mBackendId;
    lock.unlock();
    try {
        boost::filesystem::create_directories(GenerateDirname(storage, jobId));
    } catch (boost::filesystem::filesystem_error &exc) {
        SynchCout(exc.what());
    }
    for (IterIdx idx = minIterIdx; idx <= maxIterIdx; ++idx) {
        bool loaded;
        std::ifstream inStream;
        inStream.open(GenerateFilename(storage, jobId, idx).c_str());
        if (inStream.good()) {
            loaded = true;
            inStream.close();
        } else {
            loaded = false;
        }

        if (!loaded && full) {
            try {
                IterationSnapshot snp;
                RCF::ClientProgressPtr progress(new RCF::ClientProgress());
                progress->mTriggerMask = RCF::ClientProgress::Timer;
                progress->mTimerIntervalMs = 100;
                progress->mProgressCallback = boost::bind(FrontendCommunicator::OnProgress, _5);
                RcfClient<BackendIfc> client(mBackendEndpoint);
                client.getClientStub().getTransport().setMaxMessageLength(MESSAGE_SIZE);
                client.getClientStub().setMessageFilters(mCompressFlt);
                client.getClientStub().setClientProgressPtr(progress);
                client.getClientStub().setConnectTimeoutMs(10 * 1000);
                client.getClientStub().setRemoteCallTimeoutMs(60 * 1000);
                client.getClientStub().setAutoReconnect(true);
                client.getClientStub().setTries(3);
                snp = client.GetJobHistory(jobId, idx, loaded);
                if (loaded) {
                    WriteSnapshotToFile(
                        GenerateFilename(storage, jobId, idx), snp);
                }
            } catch (RCF::Exception &exc) {
                SynchCout(exc.getErrorString());
                loaded = false;
                break;
            }
        }

        if (loaded) {
            IterationSnapshotProxy proxy(storage, jobId, idx);
            history.push_back(proxy);
            emit SendHistory(proxy);
        }

        QCoreApplication::processEvents(); // Keep GUI responsive.
    }
}

void FrontendCommunicator::SetFingerprintSelector(JobId jobId,
    FingerprintSelector selector, std::string &password)
{
    try {
        bool success =
            RcfClient<BackendIfc>(mBackendEndpoint).SetFingerprintSelector(jobId, selector, password);
        if (success) {
            emit ReportStatus(tr("Fingerprint method changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change fingerprint method (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SetSimCoeffSelector(JobId jobId,
    SimCoeffSelector selector, std::string &password)
{
    try {
        bool success =
            RcfClient<BackendIfc>(mBackendEndpoint).SetSimCoeffSelector(jobId, selector, password);
        if (success) {
            emit ReportStatus(tr("Similarity coefficient method changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change similarity coefficient method (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SetDimRedSelector(JobId jobId,
    DimRedSelector selector, std::string &password)
{
    try {
        bool success =
            RcfClient<BackendIfc>(mBackendEndpoint).SetDimRedSelector(jobId, selector, password);
        if (success) {
            emit ReportStatus(tr("Visualization method changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change visualization method (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SetChemOperSelectors(JobId jobId,
        std::vector<ChemOperSelector> &selectors, std::string &password)
{
    std::vector<boost::int32_t> converted;
    converted.resize(selectors.size(), 0);
    for (size_t i = 0; i < selectors.size(); ++i) {
        converted[i] = selectors[i];
    }

    try {
        bool success =
            RcfClient<BackendIfc>(mBackendEndpoint).SetChemOperSelectors(jobId, converted, password);
        if (success) {
            emit ReportStatus(tr("Chemical operators changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change morphing operators (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SetParams(JobId jobId, MolpherParam &params,
        std::string &password)
{
    try {
        bool success =
            RcfClient<BackendIfc>(mBackendEndpoint).SetParams(jobId, params, password);
        if (success) {
            emit ReportStatus(tr("Algorithm parameters changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change algorithm parameters (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SetDecoys(JobId jobId,
    std::vector<MolpherMolecule> &decoys, std::string &password)
{
    try {
        RcfClient<BackendIfc> client(mBackendEndpoint);
        client.getClientStub().getTransport().setMaxMessageLength(MESSAGE_SIZE);
        client.getClientStub().setMessageFilters(mCompressFlt);
        bool success = client.SetDecoys(jobId, decoys, password);
        if (success) {
            emit ReportStatus(tr("Decoy set changed"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot change decoy set (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::AddPruned(JobId jobId,
    std::vector<MolpherMolecule> &pruned, std::string &password)
{
    try {
        RcfClient<BackendIfc> client(mBackendEndpoint);
        client.getClientStub().getTransport().setMaxMessageLength(MESSAGE_SIZE);
        client.getClientStub().setMessageFilters(mCompressFlt);
        bool success =client.AddPruned(jobId, pruned, password);
        if (success) {
            emit ReportStatus(tr("Pruning scheduled"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr(
                "Cannot schedule pruning (see backend console for details)"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::EnqueueNeighborhoodTask(NeighborhoodTask &task)
{
    try {
        RcfClient<BackendIfc> client(mBackendEndpoint);
        client.getClientStub().getTransport().setMaxMessageLength(MESSAGE_SIZE);
        client.getClientStub().setMessageFilters(mCompressFlt);
        client.EnqueueNeighborhoodTask(task);
        if (task.origin.smile.empty()) {
            emit ReportStatus(tr("Coordinates recalculation scheduled"), REPORT_TIMEOUT);
        } else {
            emit ReportStatus(tr("Neighborhood generation scheduled"), REPORT_TIMEOUT);
        }
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
        emit ReportStatus(tr("Cannot send request"), REPORT_TIMEOUT);
    }
}

void FrontendCommunicator::SkipNeighborhoodTask(boost::posix_time::ptime timestamp)
{
    try {
        RcfClient<BackendIfc>(mBackendEndpoint).SkipNeighborhoodTask(timestamp);
    } catch (RCF::Exception &exc) {
        SynchCout(exc.getErrorString());
    }
}

BackendDisconnectedHandler::BackendDisconnectedHandler(
    FrontendCommunicator *communicator) : mCommunicator(communicator)
{
}

void BackendDisconnectedHandler::operator()(RCF::RcfSession &session)
{
    mCommunicator->BackendDisconnected();
}

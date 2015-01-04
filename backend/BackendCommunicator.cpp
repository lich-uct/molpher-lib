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

#include "inout.h"
#include "molpher_interface.idl"
#include "BackendCommunicator.h"

#define MESSAGE_SIZE (1024 * 1024 * 1024) // 1GB

BackendCommunicator::BackendCommunicator(
    JobManager *jobManager, NeighborhoodTaskQueue *taskQueue
    ) :
    mTaskQueue(taskQueue),
    mJobManager(jobManager),
    mPublisher(RCF::TcpEndpoint("0.0.0.0", 50001)),
    mListener(RCF::TcpEndpoint("0.0.0.0", 50002)),
    mPubSvc(new RCF::PublishingService()),
    mFltSvc(new RCF::FilterService()),
    mCompressFlt(new RCF::ZlibStatelessCompressionFilter()),
    mConnectHandler(this),
    mDisconnectHandler(this)
{
    mPubSvc->setOnConnectCallback(mConnectHandler);
    mPubSvc->setOnDisconnectCallback(mDisconnectHandler);

    mFltSvc->addFilterFactory( RCF::FilterFactoryPtr(
        new RCF::ZlibStatelessCompressionFilterFactory()));

    mPublisher.getServerTransport().setMaxMessageLength(MESSAGE_SIZE);
    mPublisher.addService(mPubSvc);
    mPublisher.addService(mFltSvc);
    mPublisher.start();
    mPubSvc->beginPublish<FrontendIfc>();

    mListener.getServerTransport().setMaxMessageLength(MESSAGE_SIZE);
    mListener.bind<BackendIfc>(*this);
    mListener.addService(mFltSvc);
    mListener.start();

    RcfClient<FrontendIfc>& client = mPubSvc->publish<FrontendIfc>();
    client.getClientStub().setMessageFilters(mCompressFlt);

    SynchCout(std::string("Communicator initialized."));
}

BackendCommunicator::~BackendCommunicator()
{
    // no-op
}

void BackendCommunicator::Halt()
{
    mPubSvc->endPublish<FrontendIfc>();

    SynchCout(std::string("Communicator halted."));
}

void BackendCommunicator::InitClient(std::string &backendId)
{
    mJobManager->OnConnect(backendId);
}

boost::uint32_t BackendCommunicator::CreateJob(IterationSnapshot snp, std::string password)
{
    return mJobManager->CreateJob(snp, password);
}

void BackendCommunicator::WakeJob(boost::uint32_t jobId, std::string password)
{
    mJobManager->WakeJob(jobId, password);
}

void BackendCommunicator::SleepJob(boost::uint32_t jobId, std::string password)
{
    mJobManager->SleepJob(jobId, password);
}

void BackendCommunicator::RemoveJob(boost::uint32_t jobId, std::string password)
{
    mJobManager->RemoveJob(jobId, password);
}

void BackendCommunicator::ChangeJobOrder(
    boost::uint32_t jobId, boost::int32_t queuePosDiff, std::string password)
{
    mJobManager->ChangeJobOrder(jobId, queuePosDiff, password);
}

bool BackendCommunicator::ValidateJobPassword(boost::uint32_t jobId, std::string password)
{
    return mJobManager->ValidateJobPassword(jobId, password);
}

IterationSnapshot BackendCommunicator::GetJobHistory(boost::uint32_t jobId,
    boost::uint32_t iterIdx, bool &loaded)
{
    IterationSnapshot snp;
    loaded = mJobManager->GetJobHistory(jobId, iterIdx, snp);
    return snp;
}

bool BackendCommunicator::SetFingerprintSelector(
    boost::uint32_t jobId, boost::int32_t selector, std::string password)
{
    return mJobManager->SetFingerprintSelector(jobId, (FingerprintSelector) selector, password);
}

bool BackendCommunicator::SetSimCoeffSelector(
    boost::uint32_t jobId, boost::int32_t selector, std::string password)
{
    return mJobManager->SetSimCoeffSelector(jobId, (SimCoeffSelector) selector, password);
}

bool BackendCommunicator::SetDimRedSelector(
    boost::uint32_t jobId, boost::int32_t selector, std::string password)
{
    return mJobManager->SetDimRedSelector(jobId, (DimRedSelector) selector, password);
}

bool BackendCommunicator::SetChemOperSelectors(
    boost::uint32_t jobId, std::vector<boost::int32_t> selectors, std::string password)
{
    std::vector<ChemOperSelector> converted;
    converted.resize(selectors.size(), (ChemOperSelector) 0);
    for (size_t i = 0; i < selectors.size(); ++i) {
        converted[i] = (ChemOperSelector) selectors[i];
    }
    return mJobManager->SetChemOperSelectors(jobId, converted, password);
}

bool BackendCommunicator::SetParams(boost::uint32_t jobId, MolpherParam params, std::string password)
{
    return mJobManager->SetParams(jobId, params, password);
}

bool BackendCommunicator::SetDecoys(
    boost::uint32_t jobId, std::vector<MolpherMolecule> decoys, std::string password)
{
    return mJobManager->SetDecoys(jobId, decoys, password);
}

bool BackendCommunicator::AddPruned(
    boost::uint32_t jobId, std::vector<MolpherMolecule> pruned, std::string password)
{
    return mJobManager->AddPruned(jobId, pruned, password);
}

void BackendCommunicator::EnqueueNeighborhoodTask(NeighborhoodTask task)
{
    mTaskQueue->Push(task);
}

void BackendCommunicator::SkipNeighborhoodTask(boost::posix_time::ptime timestamp)
{
    mTaskQueue->SkipNeighborhoodTask(timestamp);
}

void BackendCommunicator::PublishJobs(JobGroup &jobs)
{
    try {
        mPubSvc->publish<FrontendIfc>().AcceptJobs(jobs);
    } catch (RCF::Exception &exc) {
        // no-op
    }
}

void BackendCommunicator::PublishIteration(IterationSnapshot &snp)
{
    try {
        mPubSvc->publish<FrontendIfc>().AcceptIteration(snp);
    } catch (RCF::Exception &exc) {
        // no-op
    }
}

void BackendCommunicator::PublishNeighborhoodTaskResult(NeighborhoodTaskResult &res)
{
    try {
        mPubSvc->publish<FrontendIfc>().AcceptNeighborhoodTaskResult(res);
    } catch (RCF::Exception &exc) {
        // no-op
    }
}

FrontendConnectedHandler::FrontendConnectedHandler(
    BackendCommunicator *communicator) : mCommunicator(communicator)
{
}

void FrontendConnectedHandler::operator()(
    RCF::RcfSession &session, const std::string &iface)
{
    const RCF::IpAddress &address =
        dynamic_cast<const RCF::IpAddress &>(session.getRemoteAddress());
    SynchCout(std::string("Frontend connected: ").append(address.string()));
}

FrontendDisconnectedHandler::FrontendDisconnectedHandler(
    BackendCommunicator *communicator) : mCommunicator(communicator)
{
}

void FrontendDisconnectedHandler::operator()(
    RCF::RcfSession &session, const std::string &iface)
{
    SynchCout(std::string("Frontend disconnected."));
}

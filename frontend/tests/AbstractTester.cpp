/*
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

#include <cassert>

#include <QtGui/QApplication>

#include "inout.h"
#include "AbstractTester.h"
#include "auxiliary/PasswordCache.h"

AbstractTester::AbstractTester(const std::string &name) :
    mScheduledAdvancer(ADV_DISPLAY_ONLINE_STATE)
{

    mName = name;
    mConnected = false;
    std::string logFile = "TestFiles/" + mName + ".log";

    mLogFile.open(logFile.c_str());
}

AbstractTester::~AbstractTester()
{
     mLogFile.close();
}

void AbstractTester::FillSnapshot(IterationSnapshot &snp)
{
    snp.source = mSource;
    snp.target = mTarget;

    std::vector<boost::int32_t> chemOperSelectors;
    chemOperSelectors.push_back(OP_ADD_ATOM);
    chemOperSelectors.push_back(OP_REMOVE_ATOM);
    chemOperSelectors.push_back(OP_ADD_BOND);
    chemOperSelectors.push_back(OP_REMOVE_BOND);
    chemOperSelectors.push_back(OP_MUTATE_ATOM);
    chemOperSelectors.push_back(OP_INTERLAY_ATOM);
    chemOperSelectors.push_back(OP_BOND_REROUTE);
    chemOperSelectors.push_back(OP_BOND_CONTRACTION);

    snp.chemOperSelectors = chemOperSelectors;

    snp.params.cntIterations = 1000;

    if (!snp.IsValid()) {
        Log("Snapshot setting is not valid!");
    }
}

void AbstractTester::PrintSnapshot(const IterationSnapshot &snp)
{
    std::stringstream s;
    s << std::endl;
    s << "Job Id: " << snp.jobId << std::endl;
    s << "Iteration: " << snp.iterIdx << std::endl;
    s << "Candidates count: " << snp.candidates.size() << std::endl;
    s << "Decoys count: " << snp.decoys.size() << std::endl;
    s << "Pruned in last iteration: " << snp.prunedDuringThisIter.size() << std::endl;
    s << "Elapsed seconds: " << snp.elapsedSeconds;
    Log(s.str());
}

void AbstractTester::PrintJobGroup(const JobGroup &group)
{
    std::stringstream s;
    JobGroup::JobQueue::const_iterator it;
    s << "Job count: " << group.mJobMap.size() << std::endl;

    s << "Live jobs: ";
    for (it = group.mLiveJobQueue.begin(); it != group.mLiveJobQueue.end(); ++it) {
        s << (*it) << " ";
    }
    s << std::endl;

    s << "Sleeping jobs: ";
    for (it = group.mSleepingJobQueue.begin(); it != group.mSleepingJobQueue.end(); ++it) {
        s << (*it) << " ";
    }
    s << std::endl;

    s << "Terminated jobs: ";
    for (it = group.mFinishedJobQueue.begin(); it != group.mFinishedJobQueue.end(); ++it) {
        s << (*it) << " ";
    }

    Log(s.str());
}

void AbstractTester::ReadTargetAndSource(std::string& fromFile)
{
    std::vector<MolpherMolecule> mols;
    std::string path(fromFile);

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];
}

void AbstractTester::Log(const std::string &line)
{
    mLogFile << line << std::endl;
    SynchCout(std::string(mName + ": " + line));
}

void AbstractTester::ScheduleAdvancer(Advancer advancer)
{
    mScheduledAdvancer = advancer;
}

void AbstractTester::DisconnectFromCommunicator()
{
    if (!mConnected) {
        return;
    }

    disconnect(&gCommunicator, SIGNAL(DisplayOnlineState()),
        this,  SLOT(DisplayOnlineState()));
    disconnect(&gCommunicator, SIGNAL(DisplayOfflineState()),
        this, SLOT(DisplayOfflineState()));
    disconnect(&gCommunicator, SIGNAL(UpdateJobs(const JobGroup &)),
        this, SLOT(UpdateJobs(const JobGroup &)));
    disconnect(&gCommunicator, SIGNAL(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)),
        this,  SLOT(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)));
    disconnect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshot &)),
        this, SLOT(VisualizeIteration(const IterationSnapshot &)));

    mConnected = false;
}

void AbstractTester::FinishTest(TestResult result)
{
    for (unsigned int i = 0; i < mCreatedJobs.size(); ++i) {
        std::string password;
        JobId jobId = mCreatedJobs[i];
        PasswordCache::ResolvePassword(jobId, password);
        gCommunicator.SleepJob(jobId, password);
        gCommunicator.RemoveJob(jobId, password);
    }
    emit NotifyTestCompleted((int)result);
}

void AbstractTester::ConnectToCommunicator()
{
    if (mConnected) {
        return;
    }

    qRegisterMetaType<JobGroup>();
    qRegisterMetaType<IterationSnapshot>();
    qRegisterMetaType<NeighborhoodTaskResult>();

    connect(&gCommunicator, SIGNAL(DisplayOnlineState()),
        this,  SLOT(DisplayOnlineState()), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(DisplayOfflineState()),
        this, SLOT(DisplayOfflineState()), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(UpdateJobs(const JobGroup &)),
        this, SLOT(UpdateJobs(const JobGroup &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)),
        this,  SLOT(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshot &)),
        this, SLOT(VisualizeIteration(const IterationSnapshot &)), Qt::QueuedConnection);

    mConnected = true;
}
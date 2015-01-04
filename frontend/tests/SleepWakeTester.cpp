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
#include "SleepWakeTester.h"
#include "auxiliary/PasswordCache.h"

SleepWakeTester::SleepWakeTester(const std::string &password) :
    AbstractTester(std::string("SleepWakeTester")),
    mPhase(CREATE_JOB)
{

    mPassword = password;

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];
}

SleepWakeTester::~SleepWakeTester()
{
    mLogFile.close();
}

void SleepWakeTester::DisplayOnlineState()
{
    Log("Tester online");
}

void SleepWakeTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void SleepWakeTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated :");
    PrintJobGroup(jobs);

    if (mScheduledAdvancer != ADV_UPDATE_JOBS) {
        return;
    }

    switch(mPhase) {
        case CREATE_JOB:
            if (VerifyJobAlive(mJobId, jobs)) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.SleepJob(mJobId, password);
                mPhase = SLEEP_JOB;
            }
            break;

        case SLEEP_JOB:
            if (VerifyJobSleeping(mJobId, jobs)) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.WakeJob(mJobId, password);
                mPhase = WAKE_JOB;
            }
            break;

        case WAKE_JOB:
            if (VerifyJobAlive(mJobId, jobs)) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.SleepJob(mJobId, password);
                mPhase = SLEEP_JOB_2;
            }
            break;

        case SLEEP_JOB_2:
            if (VerifyJobSleeping(mJobId, jobs)) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.RemoveJob(mJobId, password);
                mPhase = REMOVE_JOB;
            }
            break;

        case REMOVE_JOB:
            if (jobs.mJobMap.find(mJobId) == jobs.mJobMap.end()) {
                FinishTest(RESULT_OK);
            }
            break;
    }
}

void SleepWakeTester::VisualizeIteration(const IterationSnapshot &snp)
{
    Log("Iteration accepted:");
    PrintSnapshot(snp);
    //NextTest(ADV_VISUALIZE_ITERATION);
}

void SleepWakeTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted :");
}

void SleepWakeTester::Test()
{
    Log("SleepWakeTest");

    IterationSnapshot jobDef;
    FillSnapshot(jobDef);

    mJobId = 0xDEADBEEF;
    gCommunicator.CreateJob(jobDef, mPassword, mJobId);
    PasswordCache::CachePassword(mJobId, mPassword);

    ScheduleAdvancer(ADV_UPDATE_JOBS);

    std::stringstream s;
    if (mJobId != 0xDEADBEEF) {
        s << "Assigned job ID: " << mJobId << std::endl; // success
    } else {
        s << "Job ID not assigned." << std::endl; // fail
    }
    Log(s.str());
}
bool SleepWakeTester::VerifyJobAlive(JobId jobId, const JobGroup& jobs)
{
    JobGroup::JobQueue::const_iterator it =
            std::find(jobs.mLiveJobQueue.begin(), jobs.mLiveJobQueue.end(), jobId);

    return it != jobs.mLiveJobQueue.end();
}

bool SleepWakeTester::VerifyJobSleeping(JobId jobId, const JobGroup& jobs)
{
    JobGroup::JobQueue::const_iterator it =
            std::find(jobs.mSleepingJobQueue.begin(), jobs.mSleepingJobQueue.end(), jobId);

    return it != jobs.mSleepingJobQueue.end();
}




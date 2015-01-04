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
#include "ChangeSelectorsTester.h"
#include "dialog/ChemOperDialog.h"
#include "auxiliary/PasswordCache.h"

ChangeSelectorsTester::ChangeSelectorsTester(bool testRunningJob, const std::string &password) :
    AbstractTester(std::string("ChangeSelectorsTester")),
    mPhase(CHANGE_FINGERPRINT),
    mTestRunningJob(testRunningJob)
{

    mPassword = password;

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    mTargetChOpSelectors.push_back(OP_ADD_ATOM);

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];
}

ChangeSelectorsTester::~ChangeSelectorsTester()
{
    mLogFile.close();
}

void ChangeSelectorsTester::DisplayOnlineState()
{
    Log("Tester online");
}

void ChangeSelectorsTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void ChangeSelectorsTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated :");
    if (!mScheduledAdvancer == ADV_UPDATE_JOBS) {
        return;
    }

    if (jobs.mJobMap.find(mJobId) == jobs.mJobMap.end()) {
        return;
    }


    PrintJobGroup(jobs);

    if (mPhase == CHANGE_FINGERPRINT && jobs.mJobMap.find(mJobId) != jobs.mJobMap.end()) {
        std::string password;
        PasswordCache::ResolvePassword(mJobId, password);
        gCommunicator.SetFingerprintSelector(mJobId, mTargetFpSelector, password);
    }

    if (mTestRunningJob) {
        ScheduleAdvancer(ADV_VISUALIZE_ITERATION);
    } else {
        ScheduleAdvancer(ADV_UPDATE_JOBS);
        JobGroup::JobMap::const_iterator it = jobs.mJobMap.find(mJobId);
        if (it != jobs.mJobMap.end()) {
            VerifyJobState(it->second);
        }
    }

}

void ChangeSelectorsTester::VisualizeIteration(const IterationSnapshot &snp)
{
    Log("Iteration accepted:");

    if (!mScheduledAdvancer == ADV_VISUALIZE_ITERATION) {
        return;
    }

    PrintSnapshot(snp);


    VerifyJobState(snp);


    //NextTest(ADV_VISUALIZE_ITERATION);
}

void ChangeSelectorsTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted :");
}

void ChangeSelectorsTester::Test()
{
    Log("SleepWakeTest");

    IterationSnapshot jobDef;
    FillSnapshot(jobDef);

    mJobId = 0xDEADBEEF;
    gCommunicator.CreateJob(jobDef, mPassword, mJobId);
    PasswordCache::CachePassword(mJobId, mPassword);
    mCreatedJobs.push_back(mJobId);

    // If the test shouldn't use the running job, create second job
    if (!mTestRunningJob) {
        gCommunicator.CreateJob(jobDef, mPassword, mJobId);
        PasswordCache::CachePassword(mJobId, mPassword);
        mCreatedJobs.push_back(mJobId);
    }

    ScheduleAdvancer(ADV_UPDATE_JOBS);

    std::stringstream s;
    if (mJobId != 0xDEADBEEF) {
        s << "Assigned job ID: " << mJobId << std::endl; // success
    } else {
        s << "Job ID not assigned." << std::endl; // fail
    }
    Log(s.str());
}

bool ChangeSelectorsTester::VerifyJobState(const IterationSnapshot& snp)
{
    if (snp.jobId != mJobId) {
        return false;
    }

    bool result = false;
    std::stringstream sstream;
    sstream << "Phase Changed: ";

    switch(mPhase) {
        case CHANGE_FINGERPRINT:
            if ((int)snp.fingerprintSelector == mTargetFpSelector) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.SetSimCoeffSelector(mJobId, mTargetSimCoeffSelector, password);
                mPhase = CHANGE_SIM_COEFF;
                sstream << "CHANGE_SIM_COEFF";
                result = true;
            }
            break;
        case CHANGE_SIM_COEFF:
            if ((int)snp.simCoeffSelector == mTargetSimCoeffSelector) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.SetDimRedSelector(mJobId, mTargetDimRedSelector, password);
                mPhase = CHANGE_DIM_RED;
                sstream << "CHANGE_DIM_RED";
                result = true;
            }
            break;
        case CHANGE_DIM_RED:
            if ((int)snp.dimRedSelector == mTargetDimRedSelector) {
                std::string password;
                PasswordCache::ResolvePassword(mJobId, password);
                gCommunicator.SetChemOperSelectors(mJobId, mTargetChOpSelectors, password);
                mPhase = CHANGE_CHEM_OPERS;
                sstream << "CHANGE_CHEM_OPERS";
                result = true;
            }
            break;
        case CHANGE_CHEM_OPERS:
            if (VerifyChemOpers(snp, mTargetChOpSelectors)) {
                FinishTest(RESULT_OK);
                sstream << "RESULT_OK";
                result = true;
            }
    }

    if (result) {
        SynchCout(sstream.str());
    }
    else {
        SynchCout("PHASE NOT CHANGED");
    }

    return result;
}

bool ChangeSelectorsTester::VerifyChemOpers(const IterationSnapshot& snp,
    std::vector<ChemOperSelector>& selectors)
{
    if (snp.chemOperSelectors.size() != selectors.size()) {
        return false;
    }

    std::vector<ChemOperSelector>::iterator it;
    std::vector<boost::int32_t>::const_iterator targetIt;
    for (it = selectors.begin(); it != selectors.end(); ++it) {
        targetIt = std::find(snp.chemOperSelectors.begin(), snp.chemOperSelectors.end(), (int)(*it) );
        if (targetIt == snp.chemOperSelectors.end()) {
            return false;
        }
    }

    return true;
}

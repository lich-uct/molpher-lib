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
#include "BasicTester.h"

#define ITERATION_MAX 1000

BasicTester::BasicTester() :
    AbstractTester(std::string("CreateJobTester"))
{
    mPassword = "pass";

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];
}

BasicTester::~BasicTester()
{
    mLogFile.close();
}

void BasicTester::DisplayOnlineState()
{
    Log("Tester online");
}

void BasicTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void BasicTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated :");
    PrintJobGroup(jobs);
}

void BasicTester::VisualizeIteration(const IterationSnapshot &snp)
{
    if (ShouldAcceptSnapshot(snp)) {
        Log("Iteration accepted:");
        PrintSnapshot(snp);


        if (mScheduledAdvancer == ADV_VISUALIZE_ITERATION) {

            if (!ValidateIteration(snp)) {
                FinishTest(RESULT_ERROR);
                return;
            }

            if (snp.iterIdx >= ITERATION_MAX)
                FinishTest(RESULT_OK);

            mLastSnapshot = snp;
        }
    }
    else {
        Log("Iteration NOT accepted:");
        PrintSnapshot(snp);
    }
}

void BasicTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted :");
}

void BasicTester::Test()
{
    Log("CreateJobTest");

    IterationSnapshot jobDef;
    FillSnapshot(jobDef);

    mLastSnapshot = jobDef;

    mJobId = 0xDEADBEEF;
    gCommunicator.CreateJob(jobDef, mPassword, mJobId);
    ScheduleAdvancer(ADV_VISUALIZE_ITERATION);

    std::stringstream s;
    if (mJobId != 0xDEADBEEF) {
        s << "Assigned job ID: " << mJobId << std::endl; // success
    } else {
        s << "Job ID not assigned." << std::endl; // fail
    }
    Log(s.str());
}

bool BasicTester::ShouldAcceptSnapshot(const IterationSnapshot& snp)
{
    if (snp.jobId != mJobId) {
       return false;
    }

    return true;
}

bool BasicTester::ValidateIteration(const IterationSnapshot& snp)
{
    if (snp.jobId != mJobId) {
        std::stringstream s;
        s << "ERROR: Bad jobId: " << snp.jobId;
        Log(s.str());
        return false;
    }

    if (snp.target.smile != mTarget.smile) {
        std::stringstream s;
        s << "ERROR: Bad target: " << snp.target.smile;
        Log(s.str());
        return false;
    }

    if (snp.source.smile != mSource.smile) {
        std::stringstream s;
        s << "ERROR: Bad source: " << snp.source.smile;
        Log(s.str());
        return false;
    }

    if (snp.iterIdx > 1 && snp.prunedDuringThisIter.empty()) {
        Log("INFO: Nothing pruned");
    }

    bool containsAll = false;;
    if (snp.candidates.size() == mLastSnapshot.candidates.size()) {
        Log("WARNING: Same count of candidates as in the last iteration");
        containsAll = true;
    }

    int shared = 0;
    double closest = 1.0;
    int leafs = 0;
    IterationSnapshot::CandidateMap::const_iterator it;
    for (it = snp.candidates.begin(); it != snp.candidates.end(); ++it) {

        // Check whether the candidate was already in last iteration
        if (mLastSnapshot.candidates.find(it->first) != mLastSnapshot.candidates.end()) {
            shared++;
        }
        else {
            containsAll = false;
        }

        // Check whether the candidate exceeded itTreshold
        if (it->second.itersWithoutDistImprovement > snp.params.itThreshold + 1) {
            std::stringstream s("ERROR: Iteration treshold violated - ");
            s << "Treshold: " << snp.params.itThreshold ;
            s << " itWithoutImprovement: " << it->second.itersWithoutDistImprovement;
            Log(s.str());
        }

        if (it->second.distToTarget < closest) {
            closest = it->second.distToTarget;
        }

        if (it->second.descendants.empty()) {
            ++leafs;
        }
    }

    {
        std::stringstream s;
        s << "INFO: Closest distance to the target: " << closest;
        Log(s.str());
    }

    {
        std::stringstream s;
        s << "INFO: Candidate tree leaf count: " << leafs;
        Log(s.str());
    }

    {
        std::stringstream s;
        s << "INFO: Candidates shared with last iteration: " << shared;
        Log(s.str());
    }

    if (containsAll) {
        Log("WARNING: Same candidates set as in last iteration");
    }

    return true;
}

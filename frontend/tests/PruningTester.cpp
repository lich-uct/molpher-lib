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
#include <vector>

#include <QtGui/QApplication>

#include "inout.h"
#include "PruningTester.h"
#include "test_utils.h"
#include "auxiliary/PasswordCache.h"

#define ITERATION_TRESHOLD 3

PruningTester::PruningTester(const std::string &password) :
    AbstractTester(std::string("PruningTester")),
    mPhase(FIRST_ITERATION),
    mBadIterationCntr(0)
{
    mPassword = password;

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];
}

PruningTester::~PruningTester()
{
    mLogFile.close();
}

void PruningTester::DisplayOnlineState()
{
    Log("Tester online");
}

void PruningTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void PruningTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated");
}

void PruningTester::VisualizeIteration(const IterationSnapshot &snp)
{
    if (ShouldAcceptSnapshot(snp)) {
        Log("Iteration accepted:");

        if (mScheduledAdvancer == ADV_VISUALIZE_ITERATION) {

            std::vector<MolpherMolecule> pruneMols;
            std::string password;

            switch(mPhase) {
                case FIRST_ITERATION:
                    if (snp.candidates.empty()) {
                        Log("WARNING: Candidate map empty");
                        break;
                    }
                    mPrunedMol =  snp.candidates.begin()->second;
                    pruneMols.push_back(mPrunedMol);
                    Log("Pruned molecule: " + mPrunedMol.smile);
                    PasswordCache::ResolvePassword(mJobId, password);
                    gCommunicator.AddPruned(mJobId, pruneMols, password);
                    mPhase = WAIT_FOR_PRUNED;
                    break;

                case WAIT_FOR_PRUNED:
                    PrintSnapshot(snp);
                    IterationSnapshot::PrunedMoleculeVector::const_iterator it;
                    Log("Pruned Mols:");
                    it = snp.prunedDuringThisIter.begin();
                    while (it != snp.prunedDuringThisIter.end()) {
                        Log(*it);
                        it++;
                    }

                    it = std::find(snp.prunedDuringThisIter.begin(),
                        snp.prunedDuringThisIter.end(), mPrunedMol.smile);
                    if (it != snp.prunedDuringThisIter.end()) {
                        Log("Molecule pruned, Test OK");
                        FinishTest(RESULT_OK);
                        break;
                    } else {
                        mBadIterationCntr++;

                        if (snp.candidates.find(mPrunedMol.smile) != snp.candidates.end()) {
                            Log("WARNING: Candidate map still contains molecule marked for pruning");
                        } else {
                            Log("INFO: Candidate map no longer contains molecule marked for pruning");
                        }

                        if (mBadIterationCntr > ITERATION_TRESHOLD) {
                            FinishTest(RESULT_ERROR);
                        }
                        break;
                    }

            }

        }
    }
}

void PruningTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted");
}

void PruningTester::Test()
{
    Log("PruningTest");

    IterationSnapshot jobDef;
    FillSnapshot(jobDef);

    mJobId = 0xDEADBEEF;
    gCommunicator.CreateJob(jobDef, mPassword, mJobId);
    PasswordCache::CachePassword( mJobId, mPassword);
    mCreatedJobs.push_back(mJobId);
    ScheduleAdvancer(ADV_VISUALIZE_ITERATION);

    std::stringstream s;
    if (mJobId != 0xDEADBEEF) {
        s << "Assigned job ID: " << mJobId << std::endl; // success
    } else {
        s << "Job ID not assigned." << std::endl; // fail
    }
    Log(s.str());
}

bool PruningTester::ShouldAcceptSnapshot(const IterationSnapshot& snp)
{
    if (snp.jobId != mJobId) {
       return false;
    }

    return true;
}
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
#include <sstream>

#include "inout.h"
#include "NeighborhoodTester.h"
#include "auxiliary/PasswordCache.h"

NeighborhoodTester::NeighborhoodTester() :
    AbstractTester(std::string("NeighborhoodTester")),
    mAttemptCount(12)
{

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 1);
    mSource = mols[0];
}

NeighborhoodTester::~NeighborhoodTester()
{
    mLogFile.close();
}

void NeighborhoodTester::DisplayOnlineState()
{
    Log("Tester online");
}

void NeighborhoodTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void NeighborhoodTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated");
}

void NeighborhoodTester::VisualizeIteration(const IterationSnapshot &snp)
{
    Log("Iteration accepted:");
}

void NeighborhoodTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted :");

    if (!ShouldAcceptNeighborhood(res)) {
        return;
    }

    if (ValidateNeighborhoodResult(res)) {
        FinishTest(RESULT_OK);
    } else {
        FinishTest(RESULT_ERROR);
    }

}

void NeighborhoodTester::Test()
{
    Log("NeighborhoodTest");

    ScheduleAdvancer(ADV_VISUALIZE_NEIGHBORHOOD_TASK_RESULT);
    NeighborhoodTask task;
    task.origin = mSource;
    task.attemptCount = mAttemptCount;
    task.maxDepth = 12;
    task.chemOperSelectors.push_back(OP_ADD_ATOM);
    task.chemOperSelectors.push_back(OP_REMOVE_ATOM);

    if (!task.IsValid()) {
        SynchCout("Definition Error: Created task is invalid!");
    }
    gCommunicator.EnqueueNeighborhoodTask(task);

}

bool NeighborhoodTester::ValidateNeighborhoodResult(const NeighborhoodTaskResult& res)
{
    if (res.origin.smile != mSource.smile) {
        return false;
    }


    if (res.reducedNeighborhood.size() != mAttemptCount) {
        SynchCout("Warning: Size and attempt count differs");
        std::stringstream s;
        s << "Size: " << res.reducedNeighborhood.size() << " Attempt: " << mAttemptCount;
        SynchCout(s.str());
    }

    return true;
}

bool NeighborhoodTester::ShouldAcceptNeighborhood(const NeighborhoodTaskResult& res)
{
    if (res.origin.smile == mSource.smile) {
       return true;
    }

    return false;
}



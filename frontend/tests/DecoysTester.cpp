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
#include "DecoysTester.h"
#include "test_utils.h"
#include "auxiliary/PasswordCache.h"

#define ITERATION_MAX 5

DecoysTester::DecoysTester(const std::string &password) :
    AbstractTester(std::string("DecoysTester")),
    mPhase(PHASE_NO_DECOYS)
{
    mPassword = password;

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");
    std::string decoysPath("TestFiles/Decoys.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 2);
    mSource = mols[0];
    mTarget = mols[1];

    ReadMolphMolsFromFile(decoysPath, mDecoys);
    std::stringstream s;
    s << "Decoys count: " << mDecoys.size();
    Log(s.str());

}

DecoysTester::~DecoysTester()
{
    mLogFile.close();
    std::map<int, std::vector<double> *>::iterator it;

    for (it = mNoDecoysDistances.begin(); it != mNoDecoysDistances.end(); ++it) {
        delete it->second;
    }

    for (it = mDecoysDistances.begin(); it != mDecoysDistances.end(); ++it) {
        delete it->second;
    }

}

void DecoysTester::DisplayOnlineState()
{
    Log("Tester online");
}

void DecoysTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void DecoysTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated");
}

void DecoysTester::VisualizeIteration(const IterationSnapshot &snp)
{
    if (ShouldAcceptSnapshot(snp)) {
        Log("Iteration accepted:");

        if (mScheduledAdvancer == ADV_VISUALIZE_ITERATION) {

            Log("Computing distances");
            ComputeDistances(snp);
            Log("Distances computed");

            if (snp.iterIdx >= ITERATION_MAX) {
                switch(mPhase) {
                    case PHASE_NO_DECOYS:
                        gCommunicator.SleepJob(mJobId, mPassword);
                        CreateDecoysJob();
                        mPhase = PHASE_DECOYS;
                        break;

                    case PHASE_DECOYS:
                        gCommunicator.SleepJob(mJobId, mPassword);
                        mPhase = PHASE_COMPUTING;
                        PrintResults();
                        FinishTest(RESULT_OK);
                        break;
                    default:
                        break;
                }
            }
        }
    }
    else {
        Log("Iteration NOT accepted:");
    }
}

void DecoysTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted");
}

void DecoysTester::Test()
{
    Log("DecoysTest");

    IterationSnapshot jobDef;
    FillSnapshot(jobDef);

    mLastSnapshot = jobDef;

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

bool DecoysTester::ShouldAcceptSnapshot(const IterationSnapshot& snp)
{
    if (snp.jobId != mJobId) {
       return false;
    }

    return true;
}

void DecoysTester::PrintResults()
{
    assert(mNoDecoysDistances.size() == mDecoysDistances.size());

    int noDecoysCloser = 0;
    int decoysCloser = 0;
    int sameDistances = 0;

    for (int i = 0; i < ITERATION_MAX; ++i) {
        std::stringstream s;
        s << "-----------------" << "ITERATION: " << i << "-----------------";
        Log(s.str());

        std::vector<double> &noDecoysDist = *(mNoDecoysDistances.find(i+1)->second);
        std::vector<double> &decoysDist = *(mDecoysDistances.find(i+1)->second);

        for (unsigned int j = 0; j < decoysDist.size(); ++j) {
            std::stringstream d;
            d << "Decoy: " << j << ": ";
            if (noDecoysDist[j] > decoysDist[j]) {
                d << " WITH decoys closer";
                decoysCloser++;
            } else if ( noDecoysDist[j] == decoysDist[j]) {
                d << " SAME distances";
                sameDistances++;
            } else {
                d << " NO decoys closer";
                noDecoysCloser++;
            }

            Log(d.str());

        }
    }

    {
        std::stringstream s;
        s << "Total WITH decoys closer: " << decoysCloser;
        Log(s.str());
    }

    {
        std::stringstream s;
        s << "Total NO decoys closer: " << noDecoysCloser;
        Log(s.str());
    }

    {
        std::stringstream s;
        s << "Total sameDistances: " << sameDistances;
        Log(s.str());
    }

}

void DecoysTester::ComputeDistances(const IterationSnapshot& snp)
{
    std::vector<double> *distances = new std::vector<double>;

    for(unsigned int i = 0; i < mDecoys.size(); ++i) {
        IterationSnapshot::CandidateMap::const_iterator it = snp.candidates.begin();

        // TODO: is dist 1.0 max?
        double minDist = 1.0;

        while (it != snp.candidates.end()) {


            double dist = ComputeDistance(it->second, mDecoys[i]);
            if (dist < minDist) {
                minDist = dist;
            }

            it++;
        }

        distances->push_back(minDist);
    }

    switch(mPhase) {
        case PHASE_NO_DECOYS:
            mNoDecoysDistances.insert(std::make_pair(snp.iterIdx, distances));
            break;

        case PHASE_DECOYS:
            mDecoysDistances.insert(std::make_pair(snp.iterIdx, distances));
            break;
        default:
            break;
    }

}

void DecoysTester::CreateDecoysJob()
{
    IterationSnapshot snp;
    FillSnapshot(snp);
    snp.decoys = mDecoys;
    gCommunicator.CreateJob(snp, mPassword, mJobId);
    PasswordCache::CachePassword(mJobId, mPassword);
    mCreatedJobs.push_back(mJobId);
}
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
#include "DatabaseTester.h"

DatabaseTester::DatabaseTester() :
    AbstractTester(std::string("DatabaseTester"))
{
    mSimilarityTester = new PubchemSimilarityTester();

    std::vector<MolpherMolecule> mols;
    std::string path("TestFiles/CID_10635-CID_15951529.sdf");

    ReadMolphMolsFromFile(path, mols);
    assert(mols.size () >= 1);
    mSource = mols[0];
}

DatabaseTester::~DatabaseTester()
{
    mLogFile.close();
    delete mSimilarityTester;
}

void DatabaseTester::DisplayOnlineState()
{
    Log("Tester online");
}

void DatabaseTester::DisplayOfflineState()
{
    Log("Tester offline");
}

void DatabaseTester::UpdateJobs(const JobGroup &jobs)
{
    Log("Jobs updated");
}

void DatabaseTester::VisualizeIteration(const IterationSnapshot &snp)
{
    Log("Iteration accepted");
}

void DatabaseTester::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{
    Log("Neighborhood result accepted");
}

void DatabaseTester::Test()
{
    Log("DatabaseTest");
    std::list<std::string> smileList;
    smileList.push_back(mSource.smile);
    mSimilarityTester->SimilaritySearch(smileList);
}
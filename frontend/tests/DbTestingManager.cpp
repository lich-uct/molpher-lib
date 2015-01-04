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

#include "components/ChemicalSpaceView.h"
#include "PubchemSimilarityTester.h"
#include "DbTestingManager.h"

DbTestingManager::DbTestingManager()
{
    // no-op
}

DbTestingManager::~DbTestingManager()
{
    DeleteCompletedTesters();
}

DbTestingManager DbTestingManager::instance;

void DbTestingManager::PerformIdentitySearch(const std::list<std::string> &smileList,
    DatabaseSelector database, ChemicalSpaceView *consument)
{
    instance.DeleteCompletedTesters();
    SimilarityTester *tester = NULL;
    switch (database) {
    case DB_PUBCHEM:
        tester = new PubchemSimilarityTester();
        break;
    case DB_CHEMBL:
        // Not implemented yet
        assert(false);
        break;
    }

    QObject::connect(tester,
        SIGNAL(PublishIdentityResult(std::string, std::list<std::string> &)),
        consument,
        SLOT(OnAcceptIdentitySearchResult(std::string, std::list<std::string> &)));

    tester->IdentitySearch(smileList);
}

void DbTestingManager::PerformSimilaritySearch(const std::list<std::string> &smileList,
    DatabaseSelector database, ChemicalSpaceView *consument)
{
    instance.DeleteCompletedTesters();
    SimilarityTester *tester = NULL;
    switch (database) {
    case DB_PUBCHEM:
        tester = new PubchemSimilarityTester();
        break;
    case DB_CHEMBL:
        // Not implemented yet
        assert(false);
        break;
    }

    QObject::connect(tester,
        SIGNAL(PublishSimilarityResult(std::string, std::list<std::string> &)),
        consument,
        SLOT(OnAcceptSimilaritySearchResult(std::string, std::list<std::string> &)));

    tester->SimilaritySearch(smileList);
}

void DbTestingManager::DeleteTester(SimilarityTester *tester)
{
    instance.DeleteCompletedTesters();
    instance.mFinishedList.push_back(tester);
}

void DbTestingManager::DeleteCompletedTesters()
{
    // Delete all finished testers
    SimilarityTester *t;
    while (!mFinishedList.empty()) {
        t = mFinishedList.front();
        mFinishedList.pop_front();
        delete t;
    }
}

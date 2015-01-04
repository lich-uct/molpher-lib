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

#include "TestManager.h"
#include "inout.h"

#include <QtGui/QApplication>

#include <cassert>

TestManager::TestManager()
{

    qRegisterMetaType<JobGroup>();
    qRegisterMetaType<IterationSnapshot>();
    qRegisterMetaType<NeighborhoodTaskResult>();

    connect(&gCommunicator, SIGNAL(DisplayOnlineState(const QString &)),
        this,  SLOT(DisplayOnlineState(const QString &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(DisplayOfflineState()),
        this, SLOT(DisplayOfflineState()), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(UpdateJobs(const JobGroup &)),
        this, SLOT(UpdateJobs(const JobGroup &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)),
        this,  SLOT(VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshot &)),
        this, SLOT(VisualizeIteration(const IterationSnapshot &)), Qt::QueuedConnection);

    mCurrentTester = NULL;
}

TestManager::~TestManager()
{
    // Delete all testers in the finished list
    AbstractTester *tester;
    while(!mFinishedList.empty()) {
        tester = mFinishedList.front();
        mFinishedList.pop_front();
        delete tester;
    }
}

void TestManager::DisplayOnlineState(const QString &ip)
{
    SynchCout("TestManager: DisplayOnlineState");
    NextTest();
}

void TestManager::DisplayOfflineState()
{
    SynchCout("TestManager: DisplayOfflineState");
    QCoreApplication::exit(0); // terminate event loop, return to main()
}

void TestManager::UpdateJobs(const JobGroup &jobs)
{

}

void TestManager::VisualizeIteration(const IterationSnapshot &snp)
{

}

void TestManager::VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res)
{

}

void TestManager::AcceptTestCompleted(int result)
{
    disconnect(mCurrentTester, SIGNAL(NotifyTestCompleted(int)),
        this, SLOT(AcceptTestCompleted(int)));
    mCurrentTester->DisconnectFromCommunicator();
    mCurrentTester = NULL;

    switch ((AbstractTester::TestResult) result) {
        case AbstractTester::RESULT_OK:
            SynchCout("Test completed - Result: OK");
            break;

        case AbstractTester::RESULT_ERROR:
            SynchCout("Test completed - Result: ERROR");
            break;

        case AbstractTester::RESULT_TIMEOUT:
            SynchCout("Test completed - Result: TIMEOUT");
            break;
    }

    NextTest();
}

void TestManager::NextTest()
{
    assert(mCurrentTester == NULL);

    if (mTestList.empty()) {
        SynchCout("All tests finished");
        QCoreApplication::exit(0); // terminate event loop, return to main()
        return;
    }

    // Delete all testers in the finished list
    AbstractTester *tester;
    while(!mFinishedList.empty()) {
        tester = mFinishedList.front();
        mFinishedList.pop_front();
        delete tester;
    }

    // Move curent tester into finished list
    mFinishedList.push_back(mCurrentTester);

    // Take and run the next tester
    mCurrentTester = mTestList.front();
    mTestList.pop_front();
    mCurrentTester->ConnectToCommunicator();
    mCurrentTester->Test();
}

void TestManager::RegisterTester(AbstractTester *tester)
{
    connect(tester, SIGNAL(NotifyTestCompleted(int)),
        this, SLOT(AcceptTestCompleted(int)));

    mTestList.push_back(tester);
}

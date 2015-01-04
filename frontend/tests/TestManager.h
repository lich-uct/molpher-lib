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

#pragma once

#include "FrontendCommunicator.h"
#include "AbstractTester.h"

#include <list>

#include <QtCore/QObject>
#include <QtCore/QString>

class TestManager : public QObject
{
     Q_OBJECT
public:
    TestManager();
    ~TestManager();

    // Add tester to the test list.
    // The TestManager is responsible for connecting and disconnecting the testers.
    void RegisterTester(AbstractTester *tester);

public slots:
    void DisplayOnlineState(const QString &ip);
    void DisplayOfflineState();
    void UpdateJobs(const JobGroup &jobs);
    void VisualizeIteration(const IterationSnapshot &snp);
    void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);

    void AcceptTestCompleted(int result);

protected:
    // All testers in the list are connected to the TestManager via NotifyTestCompleted signal.
    // No tester in the list is connected to the communicator;
    std::list<AbstractTester *> mTestList;

    // List of finished testers. Tester in this list are marked "for delete"
    std::list<AbstractTester *> mFinishedList;

    // The currenly running tester. Connected to the communicator.
    // If no tester is running, this pointer is NULL
    AbstractTester *mCurrentTester;

    // Starts the next test
    void NextTest();


};


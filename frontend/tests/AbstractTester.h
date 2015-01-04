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

#include <iostream>

#include <QtCore/QObject>

#include "FrontendCommunicator.h"

class AbstractTester : public QObject
{
    Q_OBJECT
public:
    enum Advancer {
        ADV_DISPLAY_ONLINE_STATE,
        ADV_DISPLAY_OFFLINE_STATE,
        ADV_UPDATE_JOBS,
        ADV_VISUALIZE_ITERATION,
        ADV_VISUALIZE_NEIGHBORHOOD_TASK_RESULT,
        ADV_OTHER // reserved mainly for synchronous tests
    };

    enum TestResult {
        RESULT_OK,
        RESULT_TIMEOUT,
        RESULT_ERROR
    };

    void ConnectToCommunicator();
    void DisconnectFromCommunicator();
    virtual void Test()=0;
    virtual ~AbstractTester();

signals:
    void NotifyTestCompleted(int result);

protected:
    AbstractTester(const std::string &name);

    Advancer mScheduledAdvancer;
    void ScheduleAdvancer(Advancer advancer);
    void Log(const std::string &line);
    void FinishTest(TestResult result);

    void FillSnapshot(IterationSnapshot &snp);
    void PrintSnapshot(const IterationSnapshot &snp);
    void PrintJobGroup(const JobGroup &group);
    void ReadTargetAndSource(std::string &fromFile);

    std::ofstream mLogFile;

    // test variables
    std::string mPassword;
    MolpherMolecule mSource;
    MolpherMolecule mTarget;
    unsigned int mTimeout;
    std::string mName;
    bool mConnected;
    std::vector<JobId> mCreatedJobs;

public slots:
    virtual void DisplayOnlineState(){};
    virtual void DisplayOfflineState(){};
    virtual void UpdateJobs(const JobGroup &jobs){};
    virtual void VisualizeIteration(const IterationSnapshot &snp){};
    virtual void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res){};


};




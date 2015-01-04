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

#include "AbstractTester.h"

class SleepWakeTester : public AbstractTester
{
    Q_OBJECT
public:
    SleepWakeTester(const std::string &password);
    ~SleepWakeTester();

    virtual void Test();

public slots:
    virtual void DisplayOnlineState();
    virtual void DisplayOfflineState();
    virtual void UpdateJobs(const JobGroup &jobs);
    virtual void VisualizeIteration(const IterationSnapshot &snp);
    virtual void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);

protected:
    enum Phase {
        CREATE_JOB,
        SLEEP_JOB,
        WAKE_JOB,
        SLEEP_JOB_2,
        REMOVE_JOB,
    };

    JobId mJobId;
    Phase mPhase;

    bool VerifyJobAlive(JobId jobId, const JobGroup &jobs);
    bool VerifyJobSleeping(JobId jobId, const JobGroup &jobs);


};



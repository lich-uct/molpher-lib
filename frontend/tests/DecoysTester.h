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

#include "AbstractTester.h"
#include "FrontendCommunicator.h"

class DecoysTester : public AbstractTester
{
    Q_OBJECT
public:
    DecoysTester(const std::string &password);
    ~DecoysTester();

    virtual void Test();

public slots:
    virtual void DisplayOnlineState();
    virtual void DisplayOfflineState();
    virtual void UpdateJobs(const JobGroup &jobs);
    virtual void VisualizeIteration(const IterationSnapshot &snp);
    virtual void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);

private:
    enum Phase {
        PHASE_NO_DECOYS,
        PHASE_DECOYS,
        PHASE_COMPUTING
    };

    IterationSnapshot mLastSnapshot;
    JobId mJobId;
    Phase mPhase;
    std::vector<MolpherMolecule> mDecoys;
    std::map<int, std::vector<double> *> mNoDecoysDistances;
    std::map<int, std::vector<double> *> mDecoysDistances;

    bool ShouldAcceptSnapshot(const IterationSnapshot &snp);
    void PrintResults();
    void ComputeDistances(const IterationSnapshot &snp);
    void CreateDecoysJob();

};
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

class ChangeSelectorsTester : public AbstractTester {
    Q_OBJECT
public:
    ChangeSelectorsTester(bool testRunningJob, const std::string &password);
    ~ChangeSelectorsTester();

    virtual void Test();

public slots:
    virtual void DisplayOnlineState();
    virtual void DisplayOfflineState();
    virtual void UpdateJobs(const JobGroup &jobs);
    virtual void VisualizeIteration(const IterationSnapshot &snp);
    virtual void VisualizeNeighborhoodTaskResult(const NeighborhoodTaskResult &res);

protected:
    enum Phase {
        CHANGE_FINGERPRINT,
        CHANGE_SIM_COEFF,
        CHANGE_DIM_RED,
        CHANGE_CHEM_OPERS
    };

    JobId mJobId;
    Phase mPhase;

    FingerprintSelector mTargetFpSelector;
    SimCoeffSelector mTargetSimCoeffSelector;
    DimRedSelector mTargetDimRedSelector;
    std::vector<ChemOperSelector> mTargetChOpSelectors;
    bool mTestRunningJob;

    bool VerifyJobState(const IterationSnapshot &snp);
    bool VerifyChemOpers(const IterationSnapshot &snp, std::vector<ChemOperSelector> &selectors);

};



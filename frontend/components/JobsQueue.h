/*
 Copyright (c) 2012 Martin Straka
 Copyright (c) 2012 Petr Koupy

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

#include <QtGui/QDockWidget>
#include <QtGui/QTableWidget>

#include "global_types.h"
#include "meta/IterationSnapshotMeta.h"
#include "meta/JobGroupMeta.h"

class JobsQueue : public QDockWidget
{
    Q_OBJECT

public:
    JobsQueue(QWidget *parent);
    ~JobsQueue();

public slots:
    void OnUpdateJobs(const JobGroup &jobs);
    void OnVisualizeIteration(const IterationSnapshot &snp);
    void RefreshStatistics(const IterationSnapshot &snapshot);
    void CalculateStatistics(const IterationSnapshot &snapshot,
        double &distance, int &leaves);

protected slots:
    void OnShiftUp(int row);
    void OnShiftDown(int row);
    void OnAction(int row);
    void OnOpenLiveTab(int row);
    void OnOpenDetachedTab(int row);
    void OnRetrieveHistory(int row);
    void closeEvent(QCloseEvent *event);

signals:
    void ChangeJobOrder(JobId jobId, int queuePosDiff, std::string &password);
    void OpenLiveTab(JobId jobId, IterationSnapshot &snapshot);
    void OpenDetachedTab(JobId jobId, IterationSnapshot &snapshot);
    void DetachVisualizationTab(JobId jobId, IterationSnapshot &snapshot);
    void JobsUpdated(const JobGroup &jobs);
    void UpdateJobHistory(JobId jobId, IterIdx maxIterIdx);
    void CloseJobsQueue();

protected:
    enum Columns {
        COL_ID,
        COL_STATE,
        COL_SOURCE,
        COL_TARGET,
        COL_SHIFT_STEP,
        COL_SHIFT_UP,
        COL_SHIFT_DOWN,
        COL_ACTIONS,
        COL_LIVE_TAB,
        COL_DETACHED_TAB,
        COL_HISTORY,
        COL_ITERATION,
        COL_DISTANCE,
        COL_CANDIDATE_COUNT,
        COL_LEAF_COUNT,
        COL_PRUNED_COUNT,
        COL_ELAPSED_TIME,

        COL_COUNT
    };

    enum Actions {
        ACT_SET_PARAMETERS,
        ACT_SET_CHEMOPER,
        ACT_SET_DECOYS,
        ACT_SET_STATE,
        ACT_SET_PASSWORD,
        ACT_SET_FINGERPRINT,
        ACT_SET_SIMCOEFF,
        ACT_SET_DIMRED,

        ACT_COUNT
    };

private:
    QTableWidget *mTable;
    JobGroup mJobs;
    std::string mPassword;
};

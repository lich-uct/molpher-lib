/*
 Copyright (c) 2012 Marek Mikes
 Copyright (c) 2012 Martin Straka

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

#include "auxiliary/QtMocHack.h"

#include <QtGui/QAction>
#include <QtGui/QListWidget>
#include <QtGui/QWidget>
#include <QtGui/QPushButton>

#include "global_types.h"
#include "ChemicalSpaceView.h"

class Tab : public QWidget
{
    Q_OBJECT

public:
    Tab(JobId jobId, IterationSnapshot &snapshot, bool detached);
    Tab(std::vector<IterationSnapshotProxy> &snapshots);
    ~Tab();

    bool IsOnline() const;
    JobId GetJobId() const;
    void FillSnapshotsList(IterIdx minIterIdx, IterIdx maxIterIdx, bool full);
    void FillSnapshotsList(std::vector<IterationSnapshotProxy> &snapshots);

public slots:
    void OnDisplayOfflineState();
    void OnVisualizeIteration(const IterationSnapshot &snapshot);
    void OnVisualizeIteration(const IterationSnapshotProxy &snapshot);
    void OnSendHistory(const IterationSnapshotProxy &snapshot);
    void OnRevisualizeSimilarMolecules(NeighborhoodTask &task);

protected:
    void Init();
    void SetMenuOfSelectButton();
    void SetMenuOfGoToButton();
    void SetMenuOfActionButton();
    void AddSnapshotToList(const IterationSnapshotProxy &snapshot);
    QString PrintStatistics(const IterationSnapshot &snapshot);

protected slots:
    void OnSnapshotEntered(QListWidgetItem *item);
    void OnSnapshotDoubleClicked(QListWidgetItem *item);
    void OnSnapshotChanged(QListWidgetItem *item, QListWidgetItem *previous);

    // slots used by select button
    void OnButtonSelectSource();
    void OnButtonSelectTarget();
    void OnButtonSelectCandidates();
    void OnButtonSelectNewCandidates();
    void OnButtonSelectNodes();
    void OnButtonSelectLeaves();
    void OnButtonSelectDecoys();
    void OnButtonSelectNeighborhoods();
    void OnButtonSelectNeighborhoodOrigin();
    void OnButtonSelectPruned();

    // slots used by action button
    void OnButtonPrune();
    void OnButtonExport();
    void OnButtonFilter();

signals:
    void Prune(JobId jobId);
    void OpenAdhocTab(IterationSnapshotProxy &snapshot);
    void OpenAdhocTab(std::vector<IterationSnapshotProxy> &snapshots);
    void RetrieveHistory(JobId jobId, IterIdx minIterIdx, IterIdx maxIterIdx,
        bool full, std::vector<IterationSnapshotProxy> &history);
    void Select(VisualizedMolecule::ColorDefinition color);
    void Select(VisualizedMolecule::ShapeDefinition shape);
    void EnqueueNeighborhoodTask(NeighborhoodTask &task);

private:
    bool mIsOnline;
    bool mHaveJobId;
    JobId mJobId;
    bool mPruneMeansDelete;
    ChemicalSpaceView *mChemicalSpace;
    QListWidget *mSnapshotsList;
    QPushButton *mButtonSelect;
    QPushButton *mButtonGoTo;
    QPushButton *mButtonAction;
    QAction *mActionPrune;
};

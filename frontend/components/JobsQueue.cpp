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

#include <QtCore/QSignalMapper>
#include <QtCore/QStringList>
#include <QtGui/QFontMetrics>
#include <QtGui/QTableWidgetItem>
#include <QtGui/QSpinBox>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>

#include "FrontendCommunicator.h"
#include "auxiliary/PasswordCache.h"

#include "dialog/dialog_helpers.h"
#include "dialog/PasswordDialog.h"
#include "dialog/ComboBoxDialog.h"
#include "dialog/ChemOperDialog.h"
#include "dialog/ParametersDialog.h"
#include "dialog/StateDialog.h"
#include "dialog/DecoysDialog.h"

#include "JobsQueue.h"

JobsQueue::JobsQueue(QWidget *parent) :
    QDockWidget(tr("Job queue"), parent)
{
    this->setWindowTitle(tr("Job queue"));
    this->setAllowedAreas(Qt::BottomDockWidgetArea);

    mTable = NULL;
    OnUpdateJobs(mJobs);

    // signals from backend
    qRegisterMetaType<JobGroup>();
    connect(&gCommunicator, SIGNAL(UpdateJobs(const JobGroup &)),
        this, SLOT(OnUpdateJobs(const JobGroup &)), Qt::QueuedConnection);
    qRegisterMetaType<IterationSnapshot>();
    connect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshot &)),
        this, SLOT(OnVisualizeIteration(const IterationSnapshot &)), Qt::QueuedConnection);

    // signals to backend
    connect(this, SIGNAL(ChangeJobOrder(JobId, int, std::string &)),
        &gCommunicator, SLOT(ChangeJobOrder(JobId, int, std::string &)));
}

JobsQueue::~JobsQueue()
{
    delete mTable;
}

void JobsQueue::closeEvent(QCloseEvent *event)
 {
    emit CloseJobsQueue();
 }

void JobsQueue::OnUpdateJobs(const JobGroup &jobs)
{
    mJobs = jobs;
    QTableWidget *table = new QTableWidget(mJobs.mJobMap.size(), COL_COUNT);

    table->setAlternatingRowColors(true);
    table->setWindowTitle(tr("Job queue"));
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionMode(QAbstractItemView::NoSelection);

    QStringList header;
    header.insert(COL_ID, tr("Id"));
    header.insert(COL_SOURCE, tr("Source"));
    header.insert(COL_TARGET, tr("Target"));
    header.insert(COL_SHIFT_STEP, tr("Shift by"));
    header.insert(COL_SHIFT_UP, tr("Shift"));
    header.insert(COL_SHIFT_DOWN, tr("Shift"));
    header.insert(COL_ACTIONS, tr("Actions"));
    header.insert(COL_LIVE_TAB, tr("Visualization"));
    header.insert(COL_DETACHED_TAB, tr("Visualization"));
    header.insert(COL_HISTORY, tr("History"));
    header.insert(COL_STATE, tr("State"));
    header.insert(COL_ITERATION, tr("Iteration"));
    header.insert(COL_DISTANCE, tr("Distance"));
    header.insert(COL_CANDIDATE_COUNT, tr("Candidates"));
    header.insert(COL_LEAF_COUNT, tr("Leaves"));
    header.insert(COL_PRUNED_COUNT, tr("Pruned"));
    header.insert(COL_ELAPSED_TIME, tr("Time"));
    table->setHorizontalHeaderLabels(header);

    QSignalMapper *mapperShiftUp = new QSignalMapper(table);
    connect(mapperShiftUp, SIGNAL(mapped(int)), this, SLOT(OnShiftUp(int)));
    QSignalMapper *mapperShiftDown = new QSignalMapper(table);
    connect(mapperShiftDown, SIGNAL(mapped(int)), this, SLOT(OnShiftDown(int)));
    QSignalMapper *mapperActions = new QSignalMapper(table);
    connect(mapperActions, SIGNAL(mapped(int)), this, SLOT(OnAction(int)));
    QSignalMapper *mapperLiveTab = new QSignalMapper(table);
    connect(mapperLiveTab, SIGNAL(mapped(int)), this, SLOT(OnOpenLiveTab(int)));
    QSignalMapper *mapperDetachedTab = new QSignalMapper(table);
    connect(mapperDetachedTab, SIGNAL(mapped(int)), this, SLOT(OnOpenDetachedTab(int)));
    QSignalMapper *mapperRetrieveHistory = new QSignalMapper(table);
    connect(mapperRetrieveHistory, SIGNAL(mapped(int)), this, SLOT(OnRetrieveHistory(int)));

    for (unsigned int row = 0; row < mJobs.mJobMap.size(); ++row) {
        IterationSnapshot snapshot;
        mJobs.GetJob(row, snapshot);

        std::string passwordEmpty; // for convenience, try validate empty password
        bool passwordValid = false;
        if (PasswordCache::ResolvePassword(snapshot.jobId, passwordEmpty)) {
            passwordValid = true;
        } else {
            gCommunicator.ValidateJobPassword(snapshot.jobId, passwordEmpty, passwordValid);
            // valid password was also automatically added to the cache for future
        }
        if (!passwordValid) {
            // no-op
        }

        table->setItem(row, COL_ID,
            new QTableWidgetItem(QString::number(snapshot.jobId)));

        QString state;
        if (mJobs.IsLive(snapshot.jobId)) {
            if (row == 0) {
                state = tr("Running");
            } else {
                state = tr("Live");
            }
        } else if (mJobs.IsSleeping(snapshot.jobId)) {
            state = tr("Sleeping");
        } else if (mJobs.IsFinished(snapshot.jobId)) {
            state = tr("Finished");
        } else {
            assert(false);
        }
        table->setItem(row, COL_STATE, new QTableWidgetItem(state));

        table->setItem(row, COL_SOURCE,
            new QTableWidgetItem(QString::fromStdString(snapshot.source.formula)));

        table->setItem(row, COL_TARGET,
            new QTableWidgetItem(QString::fromStdString(snapshot.target.formula)));

        table->setItem(row, COL_ITERATION,
            new QTableWidgetItem(QString::number(snapshot.iterIdx)));

        double distance;
        int leaves;
        CalculateStatistics(snapshot, distance, leaves);
        table->setItem(row, COL_DISTANCE,
            new QTableWidgetItem(QString::number(distance)));
        table->setItem(row, COL_CANDIDATE_COUNT,
            new QTableWidgetItem(QString::number(snapshot.candidates.size())));
        table->setItem(row, COL_LEAF_COUNT,
            new QTableWidgetItem(QString::number(leaves)));
        table->setItem(row, COL_PRUNED_COUNT,
            new QTableWidgetItem(QString::number(snapshot.prunedDuringThisIter.size())));
        table->setItem(row, COL_ELAPSED_TIME,
            new QTableWidgetItem(QString::number(snapshot.elapsedSeconds)));

        QSpinBox *spinShiftStep = new QSpinBox();
        spinShiftStep->setValue(1);
        spinShiftStep->setMaximum(100);
        spinShiftStep->setMinimum(-100);
        table->setCellWidget(row, COL_SHIFT_STEP, spinShiftStep);

        QPushButton *buttonShiftUp = new QPushButton(tr("Up"));
        buttonShiftUp->setEnabled(mJobs.IsLive(snapshot.jobId) && row != 0);
        mapperShiftUp->setMapping(buttonShiftUp, row);
        connect(buttonShiftUp, SIGNAL(clicked()), mapperShiftUp, SLOT(map()));
        table->setCellWidget(row, COL_SHIFT_UP, buttonShiftUp);

        QPushButton *buttonShiftDown = new QPushButton(tr("Down"));
        buttonShiftDown->setEnabled(mJobs.IsLive(snapshot.jobId) && row != 0);
        mapperShiftDown->setMapping(buttonShiftDown, row);
        connect(buttonShiftDown, SIGNAL(clicked()), mapperShiftDown, SLOT(map()));
        table->setCellWidget(row, COL_SHIFT_DOWN, buttonShiftDown);

        QStringList actions;
        actions.insert(ACT_SET_PARAMETERS, tr("Set parameters"));
        actions.insert(ACT_SET_CHEMOPER, tr("Choose chem. operators"));
        actions.insert(ACT_SET_DECOYS, tr("Choose decoys"));
        actions.insert(ACT_SET_STATE, tr("Change job state"));
        actions.insert(ACT_SET_PASSWORD, tr("Elevate privileges"));
        actions.insert(ACT_SET_FINGERPRINT, tr("Set fingerprint method"));
        actions.insert(ACT_SET_SIMCOEFF, tr("Set sim. coeff. method"));
        actions.insert(ACT_SET_DIMRED, tr("Set dim. red. method"));
        QComboBox *comboActions = new QComboBox();
        comboActions->addItems(actions);
        mapperActions->setMapping(comboActions, row);
        connect(comboActions, SIGNAL(activated(int)), mapperActions, SLOT(map()));
        table->setCellWidget(row, COL_ACTIONS, comboActions);

        QPushButton *buttonLiveTab = new QPushButton(tr("Live"));
        mapperLiveTab->setMapping(buttonLiveTab, row);
        connect(buttonLiveTab, SIGNAL(clicked()), mapperLiveTab, SLOT(map()));
        table->setCellWidget(row, COL_LIVE_TAB, buttonLiveTab);

        QPushButton *buttonDetachedTab = new QPushButton(tr("Detach"));
        mapperDetachedTab->setMapping(buttonDetachedTab, row);
        connect(buttonDetachedTab, SIGNAL(clicked()), mapperDetachedTab, SLOT(map()));
        table->setCellWidget(row, COL_DETACHED_TAB, buttonDetachedTab);

        QPushButton *buttonRetrieveHistory = new QPushButton(tr("Retrieve"));
        mapperRetrieveHistory->setMapping(buttonRetrieveHistory, row);
        connect(buttonRetrieveHistory, SIGNAL(clicked()), mapperRetrieveHistory, SLOT(map()));
        table->setCellWidget(row, COL_HISTORY, buttonRetrieveHistory);
    }

    QFontMetrics font(QFont("times", 12, QFont::Bold));
    table->resizeColumnToContents(COL_ID);
    table->resizeColumnToContents(COL_STATE);
    table->setColumnWidth(COL_SOURCE, font.width("_______________"));
    table->setColumnWidth(COL_TARGET, font.width("_______________"));
    table->resizeColumnToContents(COL_ITERATION);
    table->setColumnWidth(COL_SHIFT_STEP, font.width("Shift by"));
    table->setColumnWidth(COL_SHIFT_UP, font.width("Shift"));
    table->setColumnWidth(COL_SHIFT_DOWN, font.width("Down"));
    table->resizeColumnToContents(COL_ACTIONS);
    table->setColumnWidth(COL_LIVE_TAB, font.width("Visualization"));
    table->setColumnWidth(COL_DETACHED_TAB, font.width("Visualization"));
    table->setColumnWidth(COL_HISTORY, font.width("Retrieve"));
    table->resizeColumnToContents(COL_DISTANCE);
    table->resizeColumnToContents(COL_CANDIDATE_COUNT);
    table->resizeColumnToContents(COL_LEAF_COUNT);
    table->resizeColumnToContents(COL_PRUNED_COUNT);
    table->resizeColumnToContents(COL_ELAPSED_TIME);

    setWidget(table);
    delete mTable;
    mTable = table;

    emit JobsUpdated(jobs);
}

void JobsQueue::OnShiftUp(int row)
{
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);
    QSpinBox *spinBox = qobject_cast<QSpinBox *>(mTable->cellWidget(row, COL_SHIFT_STEP));
    int move = spinBox->value();
    if (move != 0) {
        if (PasswordCache::ResolvePassword(snapshot.jobId, mPassword)) {
            emit ChangeJobOrder(snapshot.jobId, -move, mPassword);
        }
    }
}

void JobsQueue::OnShiftDown(int row)
{
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);
    QSpinBox *spinBox = qobject_cast<QSpinBox *>(mTable->cellWidget(row, COL_SHIFT_STEP));
    int move = spinBox->value();
    if (move != 0) {
        if (PasswordCache::ResolvePassword(snapshot.jobId, mPassword)) {
            emit ChangeJobOrder(snapshot.jobId, move, mPassword);
        }
    }
}

void JobsQueue::OnAction(int row)
{
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);

    QComboBox *comboBox = qobject_cast<QComboBox *>(mTable->cellWidget(row, COL_ACTIONS));
    switch (comboBox->currentIndex()) {
    case ACT_SET_CHEMOPER: {
        ChemOperDialog *chemOperDialog = new ChemOperDialog(
            snapshot.jobId, snapshot.chemOperSelectors);
        chemOperDialog->show();
        break;
    }
    case ACT_SET_DECOYS: {
        DecoysDialog *decoysDialog = new DecoysDialog(
            snapshot.jobId, snapshot.decoys);
        decoysDialog->show();
        break;
    }
    case ACT_SET_PARAMETERS: {
        ParametersDialog *paramsDialog = new ParametersDialog(
            snapshot.jobId, snapshot.params);
        paramsDialog->show();
        break;
    }
    case ACT_SET_STATE: {
        StateDialog *stateDialog = NULL;
        if (mJobs.IsLive(snapshot.jobId)) {
            if (row == 0) {
                stateDialog = new StateDialog(
                    snapshot.jobId, StateDialog::SD_SLEEP);
            } else {
                stateDialog = new StateDialog(
                    snapshot.jobId, StateDialog::SD_SLEEP | StateDialog::SD_REMOVE);
            }
        } else if (mJobs.IsSleeping(snapshot.jobId)) {
            stateDialog = new StateDialog(
                snapshot.jobId, StateDialog::SD_WAKE | StateDialog::SD_REMOVE);
        } else if (mJobs.IsFinished(snapshot.jobId)) {
            stateDialog = new StateDialog(
                snapshot.jobId, StateDialog::SD_REMOVE);
        }
        stateDialog->show();
        break;
    }
    case ACT_SET_PASSWORD: {
        PasswordDialog *passwordDialog = new PasswordDialog(snapshot.jobId);
        passwordDialog->show();
        break;
    }
    case ACT_SET_FINGERPRINT: {
        ComboBoxDialog *fingerprintDialog = new ComboBoxDialog(
            ComboBoxDialog::CD_FINGERPRINT, snapshot.jobId, snapshot.fingerprintSelector);
        fingerprintDialog->show();
        break;
    }
    case ACT_SET_SIMCOEFF: {
        ComboBoxDialog *simCoeffDialog = new ComboBoxDialog(
            ComboBoxDialog::CD_SIMCOEFF, snapshot.jobId, snapshot.simCoeffSelector);
        simCoeffDialog->show();
        break;
    }
    case ACT_SET_DIMRED: {
        ComboBoxDialog *dimRedDialog = new ComboBoxDialog(
            ComboBoxDialog::CD_DIMRED, snapshot.jobId, snapshot.dimRedSelector);
        dimRedDialog->show();
        break;
    }
    default:
        assert(false);
    }
}

void JobsQueue::OnOpenLiveTab(int row)
{
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);
    emit OpenLiveTab(snapshot.jobId, snapshot);
}

void JobsQueue::OnOpenDetachedTab(int row)
{
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);
    emit OpenDetachedTab(snapshot.jobId, snapshot);
}

void JobsQueue::OnRetrieveHistory(int row)
{
    OnOpenLiveTab(row);
    IterationSnapshot snapshot;
    mJobs.GetJob(row, snapshot);
    emit UpdateJobHistory(snapshot.jobId, snapshot.iterIdx);
}

void JobsQueue::OnVisualizeIteration(const IterationSnapshot &snapshot)
{
    if (mJobs.mJobMap.find(snapshot.jobId) != mJobs.mJobMap.end()) {
        mJobs.mJobMap[snapshot.jobId] = snapshot;
        RefreshStatistics(snapshot);
    }
}

void JobsQueue::RefreshStatistics(const IterationSnapshot &snapshot)
{
    int row = mJobs.GetIndex(snapshot.jobId);
    mTable->item(row, COL_ITERATION)->setText(
        QString::number(snapshot.iterIdx));
    mTable->item(row, COL_CANDIDATE_COUNT)->setText(
        QString::number(snapshot.candidates.size()));
    mTable->item(row, COL_PRUNED_COUNT)->setText(
        QString::number(snapshot.prunedDuringThisIter.size()));
    mTable->item(row, COL_ELAPSED_TIME)->setText(
        QString::number(snapshot.elapsedSeconds));

    double distance;
    int leaves;
    CalculateStatistics(snapshot, distance, leaves);
    mTable->item(row, COL_DISTANCE)->setText(QString::number(distance));
    mTable->item(row, COL_LEAF_COUNT)->setText(QString::number(leaves));
}

void JobsQueue::CalculateStatistics(const IterationSnapshot &snp,
    double &distance, int &leaves)
{
    distance = 1.0;
    leaves = 0;
    IterationSnapshot::CandidateMap::const_iterator it;
    for (it = snp.candidates.begin(); it != snp.candidates.end(); ++it) {

        if (it->second.distToTarget < distance) {
            distance = it->second.distToTarget;
        }

        if (it->second.descendants.empty()) {
            ++leaves;
        }
    }
}

/*
 Copyright (c) 2012 Marek Mikes
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

#include <QtCore/QList>
#include <QtGui/QAction>
#include <QtGui/QHBoxLayout>
#include <QtGui/QLabel>
#include <QtGui/QListWidgetItem>
#include <QtGui/QMenu>
#include <QtGui/QFontMetrics>

#include "inout.h"
#include "dialog/dialog_helpers.h"
#include "auxiliary/PasswordCache.h"
#include "FrontendCommunicator.h"
#include "dialog/FilterDialog.h"
#include "Tab.h"

Tab::Tab(JobId jobId, IterationSnapshot &snapshot, bool detached) :
    mIsOnline(!detached),
    mHaveJobId(true),
    mJobId(jobId),
    mPruneMeansDelete(true)
{
    mChemicalSpace = new ChemicalSpaceView(snapshot, mPruneMeansDelete);
    mSnapshotsList = new QListWidget(this);
    Init();
    if (snapshot.iterIdx > 0) {
        FillSnapshotsList(snapshot.iterIdx - 1, snapshot.iterIdx, true);
    }
    FillSnapshotsList(0, snapshot.iterIdx, false);

    if (mSnapshotsList->count() >= 2) {
        QListWidgetItem *item = mSnapshotsList->item(1);
        IterationSnapshotProxy proxy =
            item->data(Qt::UserRole).value<IterationSnapshotProxy>();
        IterationSnapshot snp = Materialize(proxy);
        mChemicalSpace->Redraw(snp);
        mChemicalSpace->SetNewSnapshot(snapshot, VisualizedMolecule::SD_NEW_CANDIDATE);
    }

    if (detached) {
        OnDisplayOfflineState();
    }
}

Tab::Tab(std::vector<IterationSnapshotProxy> &snapshots) :
    mIsOnline(false),
    mHaveJobId(false),
    mJobId(-1),
    mPruneMeansDelete(true)
{
    assert(!snapshots.empty());
    IterationSnapshot last = Materialize(snapshots.back());
    mChemicalSpace = new ChemicalSpaceView(last, mPruneMeansDelete);
    mSnapshotsList = new QListWidget(this);
    Init();
    FillSnapshotsList(snapshots);
    OnDisplayOfflineState();
}

Tab::~Tab()
{
    delete mChemicalSpace;
}

void Tab::Init()
{
    QFontMetrics font(QFont("times", 12, QFont::Bold));
    mSnapshotsList->setMouseTracking(true);
    mSnapshotsList->setMinimumWidth(font.width("Iteration 9999    "));

    connect(mSnapshotsList, SIGNAL(itemEntered(QListWidgetItem *)),
        this, SLOT(OnSnapshotEntered(QListWidgetItem *)));
    connect(mSnapshotsList, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
        this, SLOT(OnSnapshotDoubleClicked(QListWidgetItem *)));
    connect(mSnapshotsList, SIGNAL(currentItemChanged(QListWidgetItem *, QListWidgetItem *)),
        this, SLOT(OnSnapshotChanged(QListWidgetItem *, QListWidgetItem *)));

    QHBoxLayout *mainLayout = new QHBoxLayout(this);
    mainLayout->addWidget(mChemicalSpace, 100);

    QVBoxLayout *rightLayout = new QVBoxLayout();
    rightLayout->addWidget(mSnapshotsList);
    mButtonSelect = new QPushButton(tr("Select"), this);
    SetMenuOfSelectButton();
    rightLayout->addWidget(mButtonSelect);
    mButtonGoTo = new QPushButton(tr("Go to"), this);
    SetMenuOfGoToButton();
    rightLayout->addWidget(mButtonGoTo);
    mButtonAction = new QPushButton(tr("Action"), this);
    SetMenuOfActionButton();
    rightLayout->addWidget(mButtonAction);
    rightLayout->setMargin(0);
    rightLayout->setSpacing(0);

    mainLayout->addLayout(rightLayout, 1);
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);
    this->setLayout(mainLayout);

    connect(this, SIGNAL(Select(VisualizedMolecule::ColorDefinition)),
        mChemicalSpace, SLOT(OnSelect(VisualizedMolecule::ColorDefinition)));
    connect(this, SIGNAL(Select(VisualizedMolecule::ShapeDefinition)),
        mChemicalSpace, SLOT(OnSelect(VisualizedMolecule::ShapeDefinition)));
    connect(this, SIGNAL(Prune(JobId)), mChemicalSpace,
        SLOT(OnPrune(JobId)));
    connect(mChemicalSpace, SIGNAL(RevisualizeSimilarMolecules(NeighborhoodTask &)),
        this, SLOT(OnRevisualizeSimilarMolecules(NeighborhoodTask &)));

    connect(&gCommunicator, SIGNAL(DisplayOfflineState()),
        this, SLOT(OnDisplayOfflineState()), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshot &)),
        this, SLOT(OnVisualizeIteration(const IterationSnapshot &)),
        Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(VisualizeIteration(const IterationSnapshotProxy &)),
        this, SLOT(OnVisualizeIteration(const IterationSnapshotProxy &)),
        Qt::QueuedConnection);
    connect(this, SIGNAL(RetrieveHistory(JobId, IterIdx, IterIdx, bool, std::vector<IterationSnapshotProxy> &)),
        &gCommunicator, SLOT(GetJobHistory(JobId, IterIdx, IterIdx, bool, std::vector<IterationSnapshotProxy> &)));
    connect(&gCommunicator, SIGNAL(SendHistory(const IterationSnapshotProxy &)),
        this, SLOT(OnSendHistory(const IterationSnapshotProxy &)));
    connect(this, SIGNAL(EnqueueNeighborhoodTask(NeighborhoodTask &)),
        &gCommunicator, SLOT(EnqueueNeighborhoodTask(NeighborhoodTask &)));
}

void Tab::SetMenuOfSelectButton()
{
    QMenu *selectMenu = new QMenu(mButtonSelect);

    QAction *action = new QAction(tr("All"), selectMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnSelectAll()));
    selectMenu->addAction(action);

    action = new QAction(tr("Source"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectSource()));
    selectMenu->addAction(action);

    action = new QAction(tr("Target"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectTarget()));
    selectMenu->addAction(action);

    action = new QAction(tr("All candidates"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectCandidates()));
    selectMenu->addAction(action);

    action = new QAction(tr("New candidates"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectNewCandidates()));
    selectMenu->addAction(action);

    action = new QAction(tr("Inner tree nodes"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectNodes()));
    selectMenu->addAction(action);

    action = new QAction(tr("Tree leaves"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectLeaves()));
    selectMenu->addAction(action);

    action = new QAction(tr("Decoys"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectDecoys()));
    selectMenu->addAction(action);

    action = new QAction(tr("Neighborhood nodes"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectNeighborhoods()));
    selectMenu->addAction(action);

    action = new QAction(tr("Neighborhood origin"), selectMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectNeighborhoodOrigin()));
    selectMenu->addAction(action);

    if (!mPruneMeansDelete) {
        action = new QAction(tr("Pruned molecules"), selectMenu);
        connect(action, SIGNAL(triggered()), this, SLOT(OnButtonSelectPruned()));
        selectMenu->addAction(action);
    }

    mButtonSelect->setMenu(selectMenu);
}

void Tab::SetMenuOfGoToButton()
{
    QMenu *goToMenu = new QMenu(mButtonGoTo);

    QAction *action = new QAction(tr("Source"), goToMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnGoToSource()));
    goToMenu->addAction(action);

    action = new QAction(tr("Target"), goToMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnGoToTarget()));
    goToMenu->addAction(action);

    action = new QAction(tr("Neighborhood origin"), goToMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnGoToNeighborhoodOrigin()));
    goToMenu->addAction(action);

    mButtonGoTo->setMenu(goToMenu);
}

void Tab::SetMenuOfActionButton()
{
    QMenu *actionMenu = new QMenu(mButtonAction);

    mActionPrune = new QAction(tr("Prune molecules"), actionMenu);
    connect(mActionPrune, SIGNAL(triggered()), this, SLOT(OnButtonPrune()));
    actionMenu->addAction(mActionPrune);

    QAction *action = new QAction(tr("Generate neighborhood"), actionMenu);
    connect(action, SIGNAL(triggered()),
        mChemicalSpace, SLOT(OnGenerateNeighborhood()));
    actionMenu->addAction(action);

    action = new QAction(tr("Recalculate coordinates"), actionMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnRevisualize()));
    actionMenu->addAction(action);

    action = new QAction(tr("Filter molecules"), actionMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonFilter()));
    actionMenu->addAction(action);

    action = new QAction(tr("Export molecules"), actionMenu);
    connect(action, SIGNAL(triggered()), this, SLOT(OnButtonExport()));
    actionMenu->addAction(action);

    action = new QAction(tr("Create molecule bookmarks"), actionMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace, SLOT(OnBookmark()));
    actionMenu->addAction(action);

    QMenu *pubchemMenu = actionMenu->addMenu(tr("Pubchem"));

    action = new QAction(tr("Match molecules"), pubchemMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace,
        SLOT(OnActionIdentityPubchem()));
    pubchemMenu->addAction(action);

    action = new QAction(tr("Generate neighborhood"), pubchemMenu);
    connect(action, SIGNAL(triggered()), mChemicalSpace,
        SLOT(OnActionSimilarityPubchem()));
    pubchemMenu->addAction(action);

    mButtonAction->setMenu(actionMenu);
}

void Tab::FillSnapshotsList(IterIdx minIterIdx, IterIdx maxIterIdx, bool full)
{
    std::vector<IterationSnapshotProxy> history;
    emit RetrieveHistory(mJobId, minIterIdx, maxIterIdx, full, history);
}

void Tab::FillSnapshotsList(std::vector<IterationSnapshotProxy> &snapshots)
{
    std::vector<IterationSnapshotProxy>::iterator it;
    for (it = snapshots.begin(); it != snapshots.end(); ++it) {
        AddSnapshotToList(*it);
    }
}

void Tab::AddSnapshotToList(const IterationSnapshotProxy &snapshot)
{
    QString snapshotTitle;
    int rowToInsert;

    if (mHaveJobId) {
        rowToInsert = 0;
        for (int i = 0; i < mSnapshotsList->count(); i++) {
            QListWidgetItem *item = mSnapshotsList->item(i);
            QVariant variant = item->data(Qt::UserRole);

            bool alreadyInList = (snapshot.iterIdx ==
                variant.value<IterationSnapshotProxy>().iterIdx);
            if (alreadyInList) {
                return;
            }

            bool positionFound = (snapshot.iterIdx >
                variant.value<IterationSnapshotProxy>().iterIdx);
            if (positionFound) {
                break;
            }

            rowToInsert++;
        }
        snapshotTitle = tr("Iteration ");
        snapshotTitle.append(QString::number(snapshot.iterIdx));
    } else {
        rowToInsert = mSnapshotsList->count();
        snapshotTitle = tr("Snapshot ");
        snapshotTitle.append(QString::number(rowToInsert));
    }

    QListWidgetItem *item = new QListWidgetItem(snapshotTitle);
    QVariant variant = QVariant::fromValue(snapshot);
    item->setData(Qt::UserRole, variant);

    mSnapshotsList->insertItem(rowToInsert, item);
}

QString Tab::PrintStatistics(const IterationSnapshot &snapshot)
{
    double distance = 1.0;
    int leaves = 0;
    IterationSnapshot::CandidateMap::const_iterator it;
    for (it = snapshot.candidates.begin(); it != snapshot.candidates.end(); ++it) {
        if (it->second.distToTarget < distance) {
            distance = it->second.distToTarget;
        }
        if (it->second.descendants.empty()) {
            ++leaves;
        }
    }

    QString summary;
    summary.append(tr("Distance=")).append(QString::number(distance)).append(" | ");
    summary.append(tr("Candidates=")).append(QString::number(snapshot.candidates.size())).append(" | ");
    summary.append(tr("Leaves=")).append(QString::number(leaves)).append(" | ");
    summary.append(tr("Pruned=")).append(QString::number(snapshot.prunedDuringThisIter.size())).append(" | ");
    summary.append(tr("Time=")).append(QString::number(snapshot.elapsedSeconds));

    return summary;
}

void Tab::OnDisplayOfflineState()
{
    mIsOnline = false;
    mActionPrune->setEnabled(false);
}

void Tab::OnSnapshotEntered(QListWidgetItem *item)
{
    QVariant variant = item->data(Qt::UserRole);
    IterationSnapshotProxy proxy = variant.value<IterationSnapshotProxy>();
    IterationSnapshot snp = Materialize(proxy);
    if (item->statusTip().isEmpty()) {
        item->setStatusTip(PrintStatistics(snp));
    }
}

void Tab::OnSnapshotDoubleClicked(QListWidgetItem *item)
{
    if (mIsOnline) {
        QVariant variant = item->data(Qt::UserRole);
        IterationSnapshotProxy proxy = variant.value<IterationSnapshotProxy>();
        emit OpenAdhocTab(proxy);
    }
}

void Tab::OnSnapshotChanged(QListWidgetItem *item, QListWidgetItem *previous)
{
    if (!mIsOnline) {
        QVariant variant = item->data(Qt::UserRole);
        IterationSnapshotProxy proxy = variant.value<IterationSnapshotProxy>();
        IterationSnapshot snp = Materialize(proxy);
        bool notFirst = mSnapshotsList->row(item) < (mSnapshotsList->count() - 1);
        if (mHaveJobId && notFirst) {
            QListWidgetItem *prevItem = mSnapshotsList->item(
                mSnapshotsList->row(item) + 1);
            IterationSnapshotProxy prevProxy =
                prevItem->data(Qt::UserRole).value<IterationSnapshotProxy>();
            IterationSnapshot prevSnp = Materialize(prevProxy);
            mChemicalSpace->Redraw(prevSnp);
            mChemicalSpace->SetNewSnapshot(snp, VisualizedMolecule::SD_NEW_CANDIDATE);
        } else {
            if (!mHaveJobId) {
                // reset path highlighting mode to mouse hover event
                mChemicalSpace->OnChangeHighlightingRule("");
            }
            mChemicalSpace->Redraw(snp);
        }
    }
}

void Tab::OnButtonSelectSource()
{
    emit Select(VisualizedMolecule::CD_SOURCE);
}

void Tab::OnButtonSelectTarget()
{
    emit Select(VisualizedMolecule::CD_TARGET);
}

void Tab::OnButtonSelectCandidates()
{
    emit Select(VisualizedMolecule::CD_NODE);
    emit Select(VisualizedMolecule::CD_LEAF);
}

void Tab::OnButtonSelectNewCandidates()
{
    emit Select(VisualizedMolecule::SD_NEW_CANDIDATE);
}

void Tab::OnButtonSelectNodes()
{
    emit Select(VisualizedMolecule::CD_NODE);
}

void Tab::OnButtonSelectLeaves()
{
    emit Select(VisualizedMolecule::CD_LEAF);
}

void Tab::OnButtonSelectDecoys()
{
    emit Select(VisualizedMolecule::CD_DECOY);
}

void Tab::OnButtonSelectNeighborhoods()
{
    emit Select(VisualizedMolecule::CD_NEIGHBOR);
}

void Tab::OnButtonSelectNeighborhoodOrigin()
{
    emit Select(VisualizedMolecule::SD_NEIGHBORHOOD_ORIGIN);
}

void Tab::OnButtonSelectPruned()
{
    emit Select(VisualizedMolecule::CD_PRUNED);
}

bool Tab::IsOnline() const
{
    return mIsOnline;
}

JobId Tab::GetJobId() const
{
    return mJobId;
}

void Tab::OnButtonPrune()
{
    assert(true == mIsOnline);

    emit Prune(mJobId);
}

void Tab::OnButtonExport()
{
    std::vector<MolpherMolecule> selection;
    mChemicalSpace->GetSelectedMolecules(selection);
    if (selection.empty()) {
        ShowWarning("There are no selected molecules for export.");
    } else {
        QString filename = GetFileNameWriteSDF(this);
        WriteMolphMolsToSDF(filename.toStdString(), selection);
    }
}

void Tab::OnButtonFilter()
{
    FilterDialog *filterDialog = new FilterDialog();
    connect(filterDialog, SIGNAL(GetSelectedMolecules(std::vector<MolpherMolecule> &)),
        mChemicalSpace, SLOT(GetSelectedMolecules(std::vector<MolpherMolecule> &)));
    filterDialog->show();
}

void Tab::OnVisualizeIteration(const IterationSnapshot &snapshot)
{
    if ((mJobId == snapshot.jobId) && mIsOnline) {
        mChemicalSpace->SetNewSnapshot(snapshot, VisualizedMolecule::SD_NEW_CANDIDATE);
    }
}

void Tab::OnVisualizeIteration(const IterationSnapshotProxy &snapshot)
{
    if ((mJobId == snapshot.jobId) && mIsOnline) {
        AddSnapshotToList(snapshot);
    }
}

void Tab::OnSendHistory(const IterationSnapshotProxy &snapshot)
{
    if (mJobId == snapshot.jobId) {
        AddSnapshotToList(snapshot);
    }
}

void Tab::OnRevisualizeSimilarMolecules(NeighborhoodTask &task)
{
    // visualization of dimension will be set, then signal is sent to frontend communicator

    QListWidgetItem *item = NULL;
    assert(mSnapshotsList->count() > 0);
    if (mIsOnline) {
        item = mSnapshotsList->item(0);
    } else {
        QList<QListWidgetItem *> selectedItems = mSnapshotsList->selectedItems();
        if (selectedItems.empty()) {
            item = mSnapshotsList->item(0);
        } else {
            item = selectedItems.first();
        }
    }
    IterationSnapshotProxy itemProxy =
        item->data(Qt::UserRole).value<IterationSnapshotProxy>();
    IterationSnapshot snp = Materialize(itemProxy);

    task.dimRedSelector = snp.dimRedSelector;
    emit EnqueueNeighborhoodTask(task);
}

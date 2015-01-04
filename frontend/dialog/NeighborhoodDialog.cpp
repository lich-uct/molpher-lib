/*
 Copyright (c) 2012 Marek Mikes

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

#include <cassert>

#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include "dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "NeighborhoodDialog.h"

const int attemtCountDefault = 50;
const int maxDepthDefault = 2;
const double maxDistanceDefault = 0.2;

NeighborhoodDialog::NeighborhoodDialog(
    MolpherMolecule &origin, std::vector<MolpherMolecule> &context) :
    QDialog(GetMainWindow()),
    mOrigin(origin),
    mContext(context),
    mDisplayFormula(true)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    UpdateMolLabels();

    FillComboWithFingerprints(mUi.comboFingerprintSelector);
    FillComboWithSimCoeffs(mUi.comboSimCoeffSelector);
    FillComboWithDimRed(mUi.comboDimRedSelector);
    LoadPersistentSettings();
    mChemOperDialog = new ChemOperDialog(mChemOperSelectors, this);

    connect(mUi.buttonChemOperSelector, SIGNAL(clicked()),
        this, SLOT(OnChemOperClick()));
    connect(mChemOperDialog, SIGNAL(SetChemOper(std::vector<ChemOperSelector> &)),
        this, SLOT(OnSetChemOper(std::vector<ChemOperSelector> &)));
    connect(mUi.buttonOk, SIGNAL(clicked()), this, SLOT(OnOk()));
    connect(mUi.buttonCancel, SIGNAL(clicked()), this, SLOT(OnCancel()));
    connect(this, SIGNAL(EnqueueNeighborhoodTask(NeighborhoodTask &)),
        &gCommunicator, SLOT(EnqueueNeighborhoodTask(NeighborhoodTask &)));
}

NeighborhoodDialog::~NeighborhoodDialog()
{
    SavePersistentSettings();
}

void NeighborhoodDialog::keyPressEvent(QKeyEvent *e)
{
    switch (e->key()) {
    case Qt::Key_Alt:
        mDisplayFormula = !mDisplayFormula;
        UpdateMolLabels();
        break;
    default:
        break;
    }
    QDialog::keyPressEvent(e);
}

void NeighborhoodDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("NeighborhoodDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    mUi.comboFingerprintSelector->setCurrentIndex(
        settings.value("fingerprintSel", DEFAULT_FP).toInt());
    mUi.comboSimCoeffSelector->setCurrentIndex(
        settings.value("simCoeffSel", DEFAULT_SC).toInt());
    mUi.comboDimRedSelector->setCurrentIndex(
        settings.value("dimRedSel", DR_KAMADAKAWAI).toInt());

    QList<QVariant> chemOpers = settings.value("chemOpers").toList();
    if (chemOpers.isEmpty()) {
        FillListWithChemOperSelectors(&mChemOperSelectors);
    } else {
        for (int i = 0; i < chemOpers.size(); ++i) {
            mChemOperSelectors.push_back(chemOpers[i].toInt());
        }
    }
    mUi.spinBoxAttemptCount->setValue(
        settings.value("attemtCount", attemtCountDefault).toInt());
    mUi.spinBoxMaxDepth->setValue(
        settings.value("maxDepth", maxDepthDefault).toInt());
    mUi.doubleSpinBoxMaxDistance->setValue(
        settings.value("maxDistance", maxDistanceDefault).toDouble());

    settings.endGroup();
}

void NeighborhoodDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("NeighborhoodDialog");
    settings.setValue("geometry", geometry());

    settings.setValue("fingerprintSel", mUi.comboFingerprintSelector->currentIndex());
    settings.setValue("simCoeffSel", mUi.comboSimCoeffSelector->currentIndex());
    settings.setValue("dimRedSel", mUi.comboDimRedSelector->currentIndex());

    QStringList chemOpers;
    for (size_t i = 0; i < mChemOperSelectors.size(); ++i) {
        chemOpers.push_back(QString::number(mChemOperSelectors[i]));
    }
    settings.setValue("chemOpers", chemOpers);

    settings.setValue("attemtCount", mUi.spinBoxAttemptCount->value());
    settings.setValue("maxDepth", mUi.spinBoxMaxDepth->value());
    settings.setValue("maxDistance", mUi.doubleSpinBoxMaxDistance->value());

    settings.endGroup();
}

void NeighborhoodDialog::UpdateMolLabels()
{
    mUi.listWidgetContext->clear();

    QString originLabel;
    if (mDisplayFormula) {
        originLabel = QString::fromStdString(mOrigin.formula);
        FillListWithFormulas(mUi.listWidgetContext, mContext);
    } else {
        originLabel = QString::fromStdString(mOrigin.smile);
        FillListWithSmiles(mUi.listWidgetContext, mContext);
    }

    if (originLabel.isEmpty()) {
        originLabel = "none";
    }
    mUi.labelOriginFormula->setText(originLabel);
}

void NeighborhoodDialog::OnChemOperClick()
{
    mChemOperDialog->show();
}

void NeighborhoodDialog::OnSetChemOper(std::vector<ChemOperSelector> &chemOperSel)
{
    mChemOperSelectors.resize(chemOperSel.size(), 0);
    for (size_t i = 0; i < chemOperSel.size(); ++i) {
        mChemOperSelectors[i] = chemOperSel[i];
    }
}

void NeighborhoodDialog::OnOk()
{
    NeighborhoodTask task;

    task.taskTimestamp = boost::posix_time::microsec_clock::universal_time();

    task.fingerprintSelector = mUi.comboFingerprintSelector->currentIndex();
    assert(-1 != task.fingerprintSelector);

    task.simCoeffSelector = mUi.comboSimCoeffSelector->currentIndex();
    assert(-1 != task.simCoeffSelector);

    task.dimRedSelector = mUi.comboDimRedSelector->currentIndex();
    assert(-1 != task.dimRedSelector);

    task.chemOperSelectors = mChemOperSelectors;

    task.origin = mOrigin;

    task.context = mContext;

    task.attemptCount = mUi.spinBoxAttemptCount->value();
    assert(task.attemptCount >= 0);

    task.maxDepth = mUi.spinBoxMaxDepth->value();
    assert(task.maxDepth >= 0);

    task.maxDistance = mUi.doubleSpinBoxMaxDistance->value();
    assert(task.maxDistance >= 0);

    if (!task.IsValid()) {
        ShowWarning(tr("Some input is invalid."));
        return;
    }

    emit SaveTimestamp(task.taskTimestamp);
    emit EnqueueNeighborhoodTask(task);

    this->accept();
}

void NeighborhoodDialog::OnCancel()
{
    this->reject();
}

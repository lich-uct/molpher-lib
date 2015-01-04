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

#include <QtCore/QSettings>
#include <QtCore/QStringList>
#include <QtGui/QFileDialog>
#include <QtGui/QDesktopWidget>

#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "dimred_selectors.h"

#include "auxiliary/GlobalObjectsHolder.h"
#include "inout.h"
#include "FrontendCommunicator.h"

#include "dialog_helpers.h"
#include "CreateSnapshotDialog.h"

CreateSnapshotDialog::CreateSnapshotDialog() : QDialog(GetMainWindow())
{
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    mUi.setupUi(this);
    FillComboWithFingerprints(mUi.comboFingerprintSelector);
    FillComboWithSimCoeffs(mUi.comboSimCoeffSelector);
    FillComboWithDimRed(mUi.comboDimRedSelector);

    LoadPersistentSettings();

    mUi.comboFingerprintSelector->setCurrentIndex(mSnapshot.fingerprintSelector);
    mUi.comboSimCoeffSelector->setCurrentIndex(mSnapshot.simCoeffSelector);
    mUi.comboDimRedSelector->setCurrentIndex(mSnapshot.dimRedSelector);

    mParamsDialog = new ParametersDialog(mSnapshot.params, this);
    mDecoysDialog = new DecoysDialog(mSnapshot.decoys, this);
    mChemOperDialog = new ChemOperDialog(mSnapshot.chemOperSelectors, this);
    mDisplayFormula = true;

    connect(mUi.buttonAddSDF, SIGNAL(clicked()),
        this, SLOT(OnButtonAddSDF()));
    connect(mUi.buttonAddBookmarks, SIGNAL(clicked()),
        this, SLOT(OnButtonAddBookmarks()));

    connect(mUi.buttonChemOperSelector, SIGNAL(clicked()),
        this, SLOT(OnButtonChemOper()));
    connect(mChemOperDialog, SIGNAL(SetChemOper(std::vector<ChemOperSelector> &)),
        this, SLOT(OnSetChemOper(std::vector<ChemOperSelector> &)));

    connect(mUi.buttonDecoysSelect, SIGNAL(clicked()),
        this, SLOT(OnButtonDecoys()));
    connect(mDecoysDialog, SIGNAL(SetDecoys(std::vector<MolpherMolecule> &)),
        this, SLOT(OnSetDecoys(std::vector<MolpherMolecule> &)));

    connect(mUi.buttonParametersChange, SIGNAL(clicked()),
        this, SLOT(OnButtonParameters()));
    connect(mParamsDialog, SIGNAL(SetParams(MolpherParam &)),
        this, SLOT(OnSetParams(MolpherParam &)));

    connect(mUi.buttonFromSavedSnapshot, SIGNAL(clicked()),
        this, SLOT(OnLoadSavedSnapshot()));
    connect(mUi.buttonSaveSnapshot, SIGNAL(clicked()),
        this, SLOT(OnSaveSnapshot()));

    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));

    connect(this, SIGNAL(CommitSnapshot(IterationSnapshot &, std::string &, JobId &)),
        &gCommunicator, SLOT(CreateJob(IterationSnapshot &, std::string &, JobId &)));
}

CreateSnapshotDialog::~CreateSnapshotDialog()
{
    SavePersistentSettings();
}

void CreateSnapshotDialog::keyPressEvent(QKeyEvent *e)
{
    switch (e->key()) {
    case Qt::Key_Alt:
        OnFormulaToggled(!mDisplayFormula);
        break;
    default:
        break;
    }
    QDialog::keyPressEvent(e);
}

void CreateSnapshotDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("CreateSnapshotDialog");

    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }

    mSnapshot.fingerprintSelector =
        settings.value("fingerprintSel", DEFAULT_FP).toInt();
    mSnapshot.simCoeffSelector =
        settings.value("simCoeffSel", DEFAULT_SC).toInt();
    mSnapshot.dimRedSelector =
        settings.value("dimRedSel", DEFAULT_DR).toInt();

    QList<QVariant> chemOpers = settings.value("chemOpers").toList();
    if (chemOpers.isEmpty()) {
        FillListWithChemOperSelectors(&mSnapshot.chemOperSelectors);
    } else {
        for (int i = 0; i < chemOpers.size(); ++i) {
            mSnapshot.chemOperSelectors.push_back(chemOpers[i].toInt());
        }
    }

    mSnapshot.params.cntCandidatesToKeep =
        settings.value("cntCandidatesToKeep", mSnapshot.params.cntCandidatesToKeep).toInt();
    mSnapshot.params.cntCandidatesToKeepMax =
        settings.value("cntCandidatesToKeepMax", mSnapshot.params.cntCandidatesToKeepMax).toInt();
    mSnapshot.params.cntIterations =
        settings.value("cntIterations", mSnapshot.params.cntIterations).toInt();
    mSnapshot.params.cntMaxMorphs =
        settings.value("cntMaxMorphs", mSnapshot.params.cntMaxMorphs).toInt();
    mSnapshot.params.cntMorphs =
        settings.value("cntMorphs", mSnapshot.params.cntMorphs).toInt();
    mSnapshot.params.cntMorphsInDepth =
        settings.value("cntMorphsInDepth", mSnapshot.params.cntMorphsInDepth).toInt();
    mSnapshot.params.distToTargetDepthSwitch =
        settings.value("distToTargetDepthSwitch", mSnapshot.params.distToTargetDepthSwitch).toDouble();
    mSnapshot.params.itThreshold =
        settings.value("itThreshold", mSnapshot.params.itThreshold).toInt();
    mSnapshot.params.timeMaxSeconds =
        settings.value("timeMaxSeconds", mSnapshot.params.timeMaxSeconds).toInt();
    mSnapshot.params.minAcceptableMolecularWeight =
        settings.value("minAcceptableMolecularWeight", mSnapshot.params.minAcceptableMolecularWeight).toDouble();
    mSnapshot.params.maxAcceptableMolecularWeight =
        settings.value("maxAcceptableMolecularWeight", mSnapshot.params.maxAcceptableMolecularWeight).toDouble();
    // TODO Load from settings? (also require save in SavePersistentSettings)
    mSnapshot.params.useSyntetizedFeasibility = true;
    
    settings.endGroup();
}

void CreateSnapshotDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("CreateSnapshotDialog");

    settings.setValue("geometry", geometry());
    settings.setValue("fingerprintSel", mSnapshot.fingerprintSelector);
    settings.setValue("simCoeffSel", mSnapshot.simCoeffSelector);
    settings.setValue("dimRedSel", mSnapshot.dimRedSelector);

    QStringList chemOpers;
    for (size_t i = 0; i < mSnapshot.chemOperSelectors.size(); ++i) {
        chemOpers.push_back(QString::number(mSnapshot.chemOperSelectors[i]));
    }
    settings.setValue("chemOpers", chemOpers);

    settings.setValue("cntCandidatesToKeep", mSnapshot.params.cntCandidatesToKeep);
    settings.setValue("cntCandidatesToKeepMax", mSnapshot.params.cntCandidatesToKeepMax);
    settings.setValue("cntIterations", mSnapshot.params.cntIterations);
    settings.setValue("cntMaxMorphs", mSnapshot.params.cntMaxMorphs);
    settings.setValue("cntMorphs", mSnapshot.params.cntMorphs);
    settings.setValue("cntMorphsInDepth", mSnapshot.params.cntMorphsInDepth);
    settings.setValue("distToTargetDepthSwitch", mSnapshot.params.distToTargetDepthSwitch);
    settings.setValue("itThreshold", mSnapshot.params.itThreshold);
    settings.setValue("timeMaxSeconds", mSnapshot.params.timeMaxSeconds);
    settings.setValue("minAcceptableMolecularWeight", mSnapshot.params.minAcceptableMolecularWeight);
    settings.setValue("maxAcceptableMolecularWeight", mSnapshot.params.maxAcceptableMolecularWeight);

    settings.endGroup();
}

void CreateSnapshotDialog::OnLoadSavedSnapshot()
{
    QString fileName = GetFileNameReadSNP(this);
    if (!fileName.isNull()) {
        IterationSnapshot snapshot;
        if (ReadSnapshotFromFile(fileName.toStdString(), snapshot)) {
            FillAllFromSnapshot(snapshot);
        }
    }
}

void CreateSnapshotDialog::FillAllFromSnapshot(IterationSnapshot &snapshot)
{
    // source and target
    mMolecules.clear();
    mMolecules.push_back(snapshot.source);    
    mMolecules.push_back(snapshot.target);

//    if (mDisplayFormula) {
        // add data to all combo boxes
        FillSpinBox(mUi.comboMoleculeSource, mMolecules);
        FillSpinBox(mUi.comboMoleculeTarget, mMolecules);
        // set indexes, we know what is in mMolecules
        mUi.comboMoleculeSource->setCurrentIndex(0);
 //   } else {
        /*
        mUi.comboMoleculeSource->clear();
        mUi.comboMoleculeSource->addItem(QString(snapshot.source.smile.c_str()));
        mUi.comboMoleculeTarget->clear();
        mUi.comboMoleculeTarget->addItem(QString(snapshot.target.smile.c_str()));
         */ 
 //   }

    // fill data structures to support job rewind
    mSnapshot = snapshot;

    //fill dialogs
    mParamsDialog->FillParamsValues(mSnapshot.params);
    mDecoysDialog->FillDecoysList(mSnapshot.decoys);
    mChemOperDialog->FillSelections(mSnapshot.chemOperSelectors);

    // selectors
    mUi.comboFingerprintSelector->setCurrentIndex(snapshot.fingerprintSelector);
    mUi.comboSimCoeffSelector->setCurrentIndex(snapshot.simCoeffSelector);
    mUi.comboDimRedSelector->setCurrentIndex(snapshot.dimRedSelector);
}

/**
 * Change whatever is formula or smile used in combo boxes.
 * @param formulaChoosen
 */
void CreateSnapshotDialog::OnFormulaToggled(bool formulaChoosen)
{
    if (formulaChoosen) {
        mDisplayFormula = true;
    } else {
        mDisplayFormula = false;
    }
    // refill data .. 
    //FillSpinBox(mUi.comboMoleculeSource, mMoleculesSource);
    //FillSpinBox(mUi.comboMoleculeTarget, mMoleculesTarget);
}

void CreateSnapshotDialog::OnButtonAddSDF()
{
    QString fileName = GetFileNameReadFile(this);
    if (!fileName.isNull()) {
        // Ask if delete old molecules?
    }
    ReadMolphMolsFromFile(fileName.toStdString(), mMolecules);
    // refresh all the combo boxes
    FillSpinBox(mUi.comboMoleculeSource, mMolecules);
    FillSpinBox(mUi.comboMoleculeTarget, mMolecules);
}
 
/**
 * Refresh data in given combo box, preserve index. Use only when adding 
 * new data!
 * @param combo
 * @param molecules Molecules to add.
 */
void CreateSnapshotDialog::FillSpinBox(QComboBox *combo, std::vector<MolpherMolecule> &molecules)
{
    // prepare data with content
    QStringList listFormulas;
    std::vector<MolpherMolecule>::iterator it;
    for (it = molecules.begin(); it != molecules.end(); ++it) {
        if (mDisplayFormula) {
            listFormulas.append(QString::fromStdString(it->formula));
        } else {
            listFormulas.append(QString::fromStdString(it->smile));
        }
    }
    // save current index
    int index = combo->currentIndex();
    // clear combo box and add data
    combo->clear();
    combo->addItems(listFormulas);
    // set index
    if (index >= 0) {
        combo->setCurrentIndex(index);
    } else {
        combo->setCurrentIndex(0);
    }
}

void CreateSnapshotDialog::OnButtonAddBookmarks()
{
    std::vector<MolpherMolecule> bookmarks;
    gGlobalObjectsHolder.GetBookmarks()->GetCurrent(bookmarks);
    // add to our molecules pool
    mMolecules.insert(mMolecules.end(), bookmarks.begin(), bookmarks.end());
    // refresh all the combo boxes
    FillSpinBox(mUi.comboMoleculeSource, mMolecules);
    FillSpinBox(mUi.comboMoleculeTarget, mMolecules);
}

void CreateSnapshotDialog::OnButtonChemOper()
{
    mChemOperDialog->show();
}

void CreateSnapshotDialog::OnSetChemOper(std::vector<ChemOperSelector> &chemOperSel)
{
    std::vector<boost::int32_t> converted;
    converted.resize(chemOperSel.size(), 0);
    for (size_t i = 0; i < chemOperSel.size(); ++i) {
        converted[i] = chemOperSel[i];
    }
    mSnapshot.chemOperSelectors = converted;
}

void CreateSnapshotDialog::OnButtonDecoys()
{
    mDecoysDialog->show();
}

void CreateSnapshotDialog::OnSetDecoys(std::vector<MolpherMolecule> &decoys)
{
    mSnapshot.decoys = decoys;
}

void CreateSnapshotDialog::OnButtonParameters()
{
    mParamsDialog->show();
}

void CreateSnapshotDialog::OnSetParams(MolpherParam &params)
{
    mSnapshot.params = params;
}

void CreateSnapshotDialog::OnOk()
{
    FillSnapshot();
    if (!mSnapshot.IsValid()) {
        ShowWarning(tr("Cannot create job, which is not fully specified."));
    } else {
        JobId jobId;
        std::string password = mUi.lineEditPassword->text().toStdString();
        for (int i = 0; i < mUi.spinJobInstances->value(); ++i) {
            emit CommitSnapshot(mSnapshot, password, jobId);
        }
    }
}

void CreateSnapshotDialog::FillSnapshot()
{
    mSnapshot.fingerprintSelector = mUi.comboFingerprintSelector->currentIndex();
    mSnapshot.simCoeffSelector = mUi.comboSimCoeffSelector->currentIndex();
    mSnapshot.dimRedSelector = mUi.comboDimRedSelector->currentIndex();
    if (mUi.comboMoleculeSource->currentIndex() >= 0) {
        mSnapshot.source = mMolecules[mUi.comboMoleculeSource->currentIndex()];
    }
    if (mUi.comboMoleculeTarget->currentIndex() >= 0) {
        mSnapshot.target = mMolecules[mUi.comboMoleculeTarget->currentIndex()];
    }
}

void CreateSnapshotDialog::OnSaveSnapshot()
{
    FillSnapshot();

    if (!mSnapshot.IsValid()) {
        ShowWarning(tr("Job is not fully specified."));
    } else {
        QString fileName = GetFileNameWriteSNP(this);
        if (!fileName.isNull()) {
            WriteSnapshotToFile(fileName.toStdString(), mSnapshot);
        }
    }
}

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

#include <string>
#include <vector>

#include <QtGui/QDialog>
#include <QtGui/QKeyEvent>

#include "global_types.h"
#include "MolpherMolecule.h"
#include "IterationSnapshot.h"

#include "ParametersDialog.h"
#include "ChemOperDialog.h"
#include "DecoysDialog.h"

#include "ui_CreateSnapshotDialog.h"

class CreateSnapshotDialog : public QDialog
{
    Q_OBJECT

public:
    CreateSnapshotDialog();
    ~CreateSnapshotDialog();

public slots:
    void OnSetParams(MolpherParam &params);
    void OnSetDecoys(std::vector<MolpherMolecule> &decoys);
    void OnSetChemOper(std::vector<ChemOperSelector> &chemOperSel);

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();
    void FillSpinBox(QComboBox *combo, std::vector<MolpherMolecule> &molecules);
    void FillAllFromSnapshot(IterationSnapshot &snapshot);
    void FillSnapshot();
    void keyPressEvent(QKeyEvent *e);

protected slots:
    void OnFormulaToggled(bool);
    void OnButtonAddSDF();
    void OnButtonAddBookmarks();
    void OnButtonDecoys();
    void OnButtonParameters();
    void OnButtonChemOper();
    void OnLoadSavedSnapshot();
    void OnSaveSnapshot();
    void OnOk();

signals:
    void CommitSnapshot(IterationSnapshot &snp, std::string &password, JobId &jobId);

private:
    Ui::CreateSnapshotDialog mUi;

    ParametersDialog *mParamsDialog;
    DecoysDialog *mDecoysDialog;
    ChemOperDialog *mChemOperDialog;

    IterationSnapshot mSnapshot;

    std::vector<MolpherMolecule> mMolecules;
    
    bool mDisplayFormula;
};

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

#pragma once

#include <QtGui/QDialog>
#include <QtGui/QKeyEvent>

#include "ChemOperDialog.h"
#include "meta/MolpherMoleculeMeta.h"
#include "meta/NeighborhoodTaskMeta.h"
#include "ui_NeighborhoodDialog.h"

class NeighborhoodDialog : public QDialog
{
    Q_OBJECT

public:
    NeighborhoodDialog(
        MolpherMolecule &origin, std::vector<MolpherMolecule> &context);
    ~NeighborhoodDialog();

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();
    void UpdateMolLabels();

protected slots:
    void OnChemOperClick();
    void OnSetChemOper(std::vector<ChemOperSelector> &chemOperSel);
    void OnOk();
    void OnCancel();
    void keyPressEvent(QKeyEvent *e);

signals:
    void EnqueueNeighborhoodTask(NeighborhoodTask &task);
    void SaveTimestamp(boost::posix_time::ptime timestamp);

private:
    Ui::NeighborhoodDialog mUi;

    std::vector<boost::int32_t> mChemOperSelectors;
    MolpherMolecule mOrigin;
    std::vector<MolpherMolecule> mContext;
    ChemOperDialog *mChemOperDialog;
    bool mDisplayFormula;
};

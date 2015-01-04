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

#include "meta/MolpherMoleculeMeta.h"
#include "meta/NeighborhoodTaskMeta.h"
#include "ui_RevisualizeDialog.h"

class RevisualizeDialog : public QDialog
{
    Q_OBJECT

public:
    RevisualizeDialog(std::vector<MolpherMolecule> &context);
    ~RevisualizeDialog();

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();

protected slots:
    void OnOk();
    void OnCancel();

signals:
    void EnqueueNeighborhoodTask(NeighborhoodTask &task);
    void SaveTimestamp(boost::posix_time::ptime timestamp);

private:
    Ui::RevisualizeDialog mUi;

    std::vector<MolpherMolecule> mContext;
};

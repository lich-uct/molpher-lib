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

#include "global_types.h"
#include "meta/MolpherMoleculeMeta.h"
#include "components/VisualizedMolecule.h"
#include "ui_DecoysDialog.h"

class DecoysDialog : public QDialog
{
    Q_OBJECT

public:
    DecoysDialog(std::vector<MolpherMolecule> &decoys, QWidget *parent);
    DecoysDialog(JobId jobId, std::vector<MolpherMolecule> &decoys);
    ~DecoysDialog();
    void FillDecoysList(std::vector<MolpherMolecule> &decoys);

protected:
    void keyPressEvent(QKeyEvent *e);

    void LoadPersistentSettings();
    void SavePersistentSettings();
    void BasicAdjustment();
    void LoadBookmarks();
    void AppendList(QListWidget *list, std::vector<MolpherMolecule> &molecules);
    void UpdateList(QListWidget *list);
    void Move(QListWidget *source, QListWidget *destination);
    bool IsIn(const QListWidgetItem *item, const QListWidget *list) const;
    void GetDecoysFromList(std::vector<MolpherMolecule> &decoys) const;
    void DeleteSelecteDecoys();

protected slots:
    void OnListOpenedDoubleClicked(QListWidgetItem *item);
    void OnListBookmarksDoubleClicked(QListWidgetItem *item);
    void OnButtonFromOpened();
    void OnButtonFromBookmarks();
    void OnButtonLoad();
    void OnButtonDeleteDecoy();
    void OnButtonSelectAllOpened();
    void OnButtonSelectAllDecoys();
    void OnButtonLoadBookmarks();
    void OnButtonSelectAllBookmarks();
    void OnOk();

signals:
    void SetDecoys(JobId jobId, std::vector<MolpherMolecule> &decoys, std::string &password);
    void SetDecoys(std::vector<MolpherMolecule> &decoys);

private:
    Ui::DecoysDialog mUi;
    int mJobId;
    bool mHaveJobId;
    bool mDisplayFormula;
};

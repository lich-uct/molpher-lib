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

#include <QtGui/QDialog>
#include <QtGui/QListWidget>

#include "chemoper_selectors.h"
#include "global_types.h"
#include "auxiliary/PasswordCache.h"
#include "ui_ChemOperDialog.h"

class ChemOperDialog : public QDialog
{
    Q_OBJECT

public:
    ChemOperDialog(std::vector<int> &selectedOper, QWidget *parent);
    ChemOperDialog(JobId jobId, std::vector<int> &selectedOper);
    ~ChemOperDialog();
    void FillSelections(std::vector<int> &selectedOper);

protected:
    void Move(QListWidget *source, QListWidget *destination);
    void SetListOperators();
    void SetConnects();
    void LoadPersistentSettings();
    void SavePersistentSettings();

protected slots:
    void OnListSetChemOperDoubleClicked(QListWidgetItem *item);
    void OnListSelectionDoubleClicked(QListWidgetItem *item);
    void OnButtonFromOperatorsSet();
    void OnButtonFromSelection();
    void OnButtonSelectAllSetChemOper();
    void OnButtonSelectAllSelection();
    void OnOk();

signals:
    void SetChemOper(std::vector<ChemOperSelector> &selectors);
    void SetChemOperSelectors(JobId jobId, std::vector<ChemOperSelector> &selectors, std::string &password);

private:
    Ui::ChemOperDialog mUi;
    JobId mJobId;
    bool mHaveJobId;
    QMap<QString, ChemOperSelector> mOperators;
};

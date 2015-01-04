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

#include <cassert>

#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include "FrontendCommunicator.h"
#include "dialog_helpers.h"
#include "ChemOperDialog.h"

ChemOperDialog::ChemOperDialog(std::vector<int> &selectedOper, QWidget *parent) :
    QDialog(parent),
    mHaveJobId(false)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);

    SetListOperators();
    FillSelections(selectedOper);
    SetConnects();
    LoadPersistentSettings();
}

ChemOperDialog::ChemOperDialog(JobId jobId, std::vector<int> &selectedOper) :
    QDialog(GetMainWindow()),
    mJobId(jobId),
    mHaveJobId(true)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    SetListOperators();
    FillSelections(selectedOper);
    SetConnects();
    LoadPersistentSettings();
}

ChemOperDialog::~ChemOperDialog()
{
    SavePersistentSettings();
}

void ChemOperDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ChemOperDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void ChemOperDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ChemOperDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void ChemOperDialog::SetListOperators()
{
    mOperators.insert(ChemOperLongDesc(OP_ADD_ATOM), OP_ADD_ATOM);
    mOperators.insert(ChemOperLongDesc(OP_REMOVE_ATOM), OP_REMOVE_ATOM);
    mOperators.insert(ChemOperLongDesc(OP_ADD_BOND), OP_ADD_BOND);
    mOperators.insert(ChemOperLongDesc(OP_REMOVE_BOND), OP_REMOVE_BOND);
    mOperators.insert(ChemOperLongDesc(OP_MUTATE_ATOM), OP_MUTATE_ATOM);
    mOperators.insert(ChemOperLongDesc(OP_INTERLAY_ATOM), OP_INTERLAY_ATOM);
    mOperators.insert(ChemOperLongDesc(OP_BOND_REROUTE), OP_BOND_REROUTE);
    mOperators.insert(ChemOperLongDesc(OP_BOND_CONTRACTION), OP_BOND_CONTRACTION);
}

void ChemOperDialog::FillSelections(std::vector<int> &selectedOper)
{
    mUi.listSelection->clear();
    mUi.listSetChemOper->clear();

    QMap<QString, ChemOperSelector>::iterator it;
    ChemOperSelector operID;
    for (it = mOperators.begin(); it != mOperators.end(); ++it) {
        operID = it.value();
        if ((selectedOper.size() != 0) &&
                (std::find(selectedOper.begin(), selectedOper.end(), operID) != selectedOper.end())) {
            mUi.listSelection->addItem(it.key());
        } else {
            mUi.listSetChemOper->addItem(it.key());
        }
    }
}

void ChemOperDialog::SetConnects()
{
    connect(mUi.listSetChemOper, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
        this, SLOT(OnListSetChemOperDoubleClicked(QListWidgetItem *)));
    connect(mUi.listSelection, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
        this, SLOT(OnListSelectionDoubleClicked(QListWidgetItem *)));
    connect(mUi.buttonFromOperatorsSet, SIGNAL(pressed()),
        this, SLOT(OnButtonFromOperatorsSet()));
    connect(mUi.buttonFromSelection, SIGNAL(pressed()),
        this, SLOT(OnButtonFromSelection()));
    connect(mUi.buttonSelectAllSetChemOper, SIGNAL(pressed()),
        this, SLOT(OnButtonSelectAllSetChemOper()));
    connect(mUi.buttonSelectAllSelection, SIGNAL(pressed()),
        this, SLOT(OnButtonSelectAllSelection()));
    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));

    connect(this, SIGNAL(SetChemOperSelectors(JobId, std::vector<ChemOperSelector> &, std::string &)),
        &gCommunicator, SLOT(SetChemOperSelectors(JobId, std::vector<ChemOperSelector> &, std::string &)));
}

void ChemOperDialog::Move(QListWidget *source, QListWidget *destination)
{
    assert(NULL != source);
    assert(NULL != destination);
    assert(source != destination);

    if ((source->selectedItems().count() == 0) && (source->currentRow() >= 0)) {
        destination->addItem(source->takeItem(source->currentRow()));
    } else {
        foreach (QListWidgetItem *item, source->selectedItems()) {
            destination->addItem(source->takeItem(source->row(item)));
        }
    }
}

void ChemOperDialog::OnOk()
{
    if (mUi.listSelection->count() > 0) {
        std::vector<ChemOperSelector> selectedChemOper;
        QListWidgetItem *item;
        ChemOperSelector chemOperSel;
        QString chemOperName;
        for (int i = 0; i < mUi.listSelection->count(); ++i) {
            item = mUi.listSelection->item(i);
            chemOperName = item->text();

            assert(mOperators.contains(chemOperName));

            chemOperSel = mOperators.value(chemOperName);
            selectedChemOper.push_back(chemOperSel);
        }

        if (mHaveJobId) {
            std::string password;
            if (PasswordCache::ResolvePassword(mJobId, password)) {
                emit SetChemOperSelectors(mJobId, selectedChemOper, password);
            } else {
                ShowPermissionError();
            }
        } else {
            emit SetChemOper(selectedChemOper);
        }
    } else {
        ShowWarning(tr("At least one morphing operator have to be chosen."));
    }
}

void ChemOperDialog::OnButtonFromOperatorsSet()
{
    Move(mUi.listSetChemOper, mUi.listSelection);
}

void ChemOperDialog::OnButtonFromSelection()
{
    Move(mUi.listSelection, mUi.listSetChemOper);
}

void ChemOperDialog::OnButtonSelectAllSetChemOper()
{
    SelectAll(mUi.listSetChemOper);
}

void ChemOperDialog::OnButtonSelectAllSelection()
{
    SelectAll(mUi.listSelection);
}

void ChemOperDialog::OnListSetChemOperDoubleClicked(QListWidgetItem *item)
{
    Q_UNUSED(item)
    Move(mUi.listSetChemOper, mUi.listSelection);
}

void ChemOperDialog::OnListSelectionDoubleClicked(QListWidgetItem *item)
{
    Q_UNUSED(item)
    Move(mUi.listSelection, mUi.listSetChemOper);
}

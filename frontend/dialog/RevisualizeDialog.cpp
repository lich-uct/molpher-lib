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
#include "RevisualizeDialog.h"

RevisualizeDialog::RevisualizeDialog(std::vector<MolpherMolecule> &context) :
    QDialog(GetMainWindow()),
    mContext(context)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    FillComboWithDimRed(mUi.comboDimRedSelector);
    mUi.comboDimRedSelector->setCurrentIndex(DR_KAMADAKAWAI);
    LoadPersistentSettings();

    connect(mUi.buttonOk, SIGNAL(clicked()), this, SLOT(OnOk()));
    connect(mUi.buttonCancel, SIGNAL(clicked()), this, SLOT(OnCancel()));
    connect(this, SIGNAL(EnqueueNeighborhoodTask(NeighborhoodTask &)),
        &gCommunicator, SLOT(EnqueueNeighborhoodTask(NeighborhoodTask &)));
}

RevisualizeDialog::~RevisualizeDialog()
{
    SavePersistentSettings();
}

void RevisualizeDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("RevisualizeDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    mUi.comboDimRedSelector->setCurrentIndex(
        settings.value("dimRedSel", DEFAULT_DR).toInt());
    settings.endGroup();
}

void RevisualizeDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("RevisualizeDialog");
    settings.setValue("geometry", geometry());
    settings.setValue("dimRedSel", mUi.comboDimRedSelector->currentIndex());
    settings.endGroup();
}

void RevisualizeDialog::OnOk()
{
    NeighborhoodTask task;

    task.taskTimestamp = boost::posix_time::microsec_clock::universal_time();

    task.dimRedSelector = mUi.comboDimRedSelector->currentIndex();
    assert(-1 != task.dimRedSelector);

    task.context = mContext;

    if (!task.IsValid()) {
        ShowWarning(tr("Some input is invalid."));
        return;
    }

    emit SaveTimestamp(task.taskTimestamp);
    emit EnqueueNeighborhoodTask(task);

    this->accept();
}

void RevisualizeDialog::OnCancel()
{
    this->reject();
}

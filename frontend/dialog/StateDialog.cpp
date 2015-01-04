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

#include <QtCore/QVariant>
#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include "auxiliary/PasswordCache.h"
#include "FrontendCommunicator.h"
#include "dialog_helpers.h"
#include "StateDialog.h"

StateDialog::StateDialog(JobId jobId, int flags) :
    QDialog(GetMainWindow()),
    mJobId(jobId)
{
    mUi.setupUi(this);
    setAttribute(Qt::WA_DeleteOnClose, true);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setWindowTitle(tr("Change job state"));
    mUi.labelDescription->setText(tr("Available actions:"));

    int itemIdx = 0;
    if (flags & SD_WAKE) {
        mUi.comboBox->insertItem(
            itemIdx++, tr("Wake the job up"), QVariant(SD_COMBO_WAKE));
    }
    if (flags & SD_SLEEP) {
        mUi.comboBox->insertItem(
            itemIdx++, tr("Put the job to sleep"), QVariant(SD_COMBO_SLEEP));
    }
    if (flags & SD_REMOVE) {
        mUi.comboBox->insertItem(
            itemIdx++, tr("Remove the job"), QVariant(SD_COMBO_REMOVE));
    }

    LoadPersistentSettings();

    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));
    connect(this, SIGNAL(WakeJob(JobId, std::string &)),
        &gCommunicator, SLOT(WakeJob(JobId, std::string &)));
    connect(this, SIGNAL(SleepJob(JobId, std::string &)),
        &gCommunicator, SLOT(SleepJob(JobId, std::string &)));
    connect(this, SIGNAL(RemoveJob(JobId, std::string &)),
        &gCommunicator, SLOT(RemoveJob(JobId, std::string &)));
}

StateDialog::~StateDialog()
{
    SavePersistentSettings();
}

void StateDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("StateDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void StateDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("StateDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void StateDialog::OnOk()
{
    std::string password;
    if (PasswordCache::ResolvePassword(mJobId, password))
    {
        int selection =
            mUi.comboBox->itemData(mUi.comboBox->currentIndex()).value<int>();
        switch (selection) {
            case SD_COMBO_WAKE:
                emit WakeJob(mJobId, password);
                break;
            case SD_COMBO_SLEEP:
                emit SleepJob(mJobId, password);
                break;
            case SD_COMBO_REMOVE:
                emit RemoveJob(mJobId, password);
                break;
            default:
                assert(false);
        }
    } else {
        ShowPermissionError();
    }
 }

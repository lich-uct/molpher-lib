/*
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

#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include "auxiliary/PasswordCache.h"
#include "dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "PasswordDialog.h"

PasswordDialog::PasswordDialog(JobId jobId) :
    QDialog(GetMainWindow()),
    mJobId(jobId)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);
    LoadPersistentSettings();
    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));
}

PasswordDialog::~PasswordDialog()
{
    SavePersistentSettings();
}

void PasswordDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("PasswordDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void PasswordDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("PasswordDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void PasswordDialog::OnOk()
{
    std::string password = mUi.editPassword->text().toStdString();
    bool isValid = false;
    if (PasswordCache::ResolvePassword(mJobId, password)) {
        isValid = true;
    } else {
        gCommunicator.ValidateJobPassword(mJobId, password, isValid);
        // valid password was also automatically added to the cache for future
    }
    if (isValid) {
        ShowInformation(tr("Password successfully verified."));
    } else {
        ShowWarning(tr("Password could not be verified."));
    }
}

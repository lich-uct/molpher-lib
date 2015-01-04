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

#include "FrontendCommunicator.h"
#include "dialog_helpers.h"
#include "ConnectDialog.h"

const int IP_COUNT = 10;

ConnectDialog::ConnectDialog() : QDialog(GetMainWindow())
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);
    mUi.comboIP->setEditable(true);
    LoadPersistentSettings();

    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));
    connect(this, SIGNAL(ConnectToBackend(std::string)),
        &gCommunicator, SLOT(ConnectToBackend(std::string)));
}

ConnectDialog::~ConnectDialog()
{
    SavePersistentSettings();
}

void ConnectDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ConnectDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }

    mListIP = settings.value("ipAddress").toStringList();
    if (mListIP.isEmpty()){
        mListIP.append("localhost");
    }
    mUi.comboIP->addItems(mListIP);
    settings.endGroup();
}

void ConnectDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ConnectDialog");
    settings.setValue("geometry", geometry());
    if (mListIP.count() > IP_COUNT) {
        mListIP.removeLast();
    }
    settings.setValue("ipAddress", mListIP);
    settings.endGroup();
}

void ConnectDialog::OnOk()
{
    std::string ip = mUi.comboIP->currentText().toStdString();
    mListIP.insert(0, QString::fromStdString(ip));
    emit ConnectToBackend(ip);
}

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

#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include "auxiliary/PasswordCache.h"
#include "dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "ComboBoxDialog.h"

ComboBoxDialog::ComboBoxDialog(int trait, JobId jobId, int currentIndex) :
    QDialog(GetMainWindow()),
    mJobId(jobId),
    mTrait(trait)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    switch (mTrait) {
    case CD_FINGERPRINT:
        this->setWindowTitle(tr("Set fingerprint method"));
        mUi.labelDescription->setText(tr("Fingerprint method:"));
        FillComboWithFingerprints(mUi.comboBox);
        break;
    case CD_SIMCOEFF:
        this->setWindowTitle(tr("Set similarity coefficient method"));
        mUi.labelDescription->setText(tr("Similarity coefficient method:"));
        FillComboWithSimCoeffs(mUi.comboBox);
        break;
    case CD_DIMRED:
        this->setWindowTitle(tr("Set visualization method"));
        mUi.labelDescription->setText(tr("Visualization method:"));
        FillComboWithDimRed(mUi.comboBox);
        break;
    default:
        assert(false);
    }

    mUi.comboBox->setCurrentIndex(currentIndex);
    LoadPersistentSettings();

    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));

    connect(this, SIGNAL(SetFingerprintSelector(JobId, FingerprintSelector, std::string &)),
        &gCommunicator, SLOT(SetFingerprintSelector(JobId, FingerprintSelector, std::string &)));
    connect(this, SIGNAL(SetSimCoeffSelector(JobId, SimCoeffSelector, std::string &)),
        &gCommunicator, SLOT(SetSimCoeffSelector(JobId, SimCoeffSelector, std::string &)));
    connect(this, SIGNAL(SetDimRedSelector(JobId, DimRedSelector, std::string &)),
        &gCommunicator, SLOT(SetDimRedSelector(JobId, DimRedSelector, std::string &)));
}

ComboBoxDialog::~ComboBoxDialog()
{
    SavePersistentSettings();
}

void ComboBoxDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ComboBoxDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void ComboBoxDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ComboBoxDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void ComboBoxDialog::OnOk()
{
    std::string password;
    if (PasswordCache::ResolvePassword(mJobId, password)) {
        switch (mTrait) {
        case CD_FINGERPRINT:
            emit SetFingerprintSelector(mJobId,
                (FingerprintSelector) mUi.comboBox->currentIndex(), password);
            break;
        case CD_SIMCOEFF:
            emit SetSimCoeffSelector(mJobId,
                (SimCoeffSelector) mUi.comboBox->currentIndex(), password);
            break;
        case CD_DIMRED:
            emit SetDimRedSelector(mJobId,
                (DimRedSelector) mUi.comboBox->currentIndex(), password);
            break;
        default:
            assert(false);
        }
    } else {
        ShowPermissionError();
    }
}

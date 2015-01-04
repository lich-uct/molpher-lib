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

#pragma once

#include <QtGui/QDialog>

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "dimred_selectors.h"
#include "simcoeff_selectors.h"
#include "ui_ComboBoxDialog.h"

class ComboBoxDialog : public QDialog
{
    Q_OBJECT

public:
    enum Trait {
        CD_FINGERPRINT,
        CD_SIMCOEFF,
        CD_DIMRED
    };

    ComboBoxDialog(int trait, JobId jobId, int currentIndex = 0);
    ~ComboBoxDialog();

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();

protected slots:
    void OnOk();

signals:
    void SetFingerprintSelector(JobId jobId, FingerprintSelector selector, std::string &password);
    void SetSimCoeffSelector(JobId jobId, SimCoeffSelector selector, std::string &password);
    void SetDimRedSelector(JobId jobId, DimRedSelector selector, std::string &password);

private:
    Ui::ComboBoxDialog mUi;
    JobId mJobId;
    int mTrait;
};

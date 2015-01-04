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

#pragma once

#include <QtGui/QDialog>

#include "global_types.h"
#include "ui_StateDialog.h"

namespace Ui {
class StateDialog;
}

class StateDialog : public QDialog
{
    Q_OBJECT

public:
     enum Flags {
         SD_EMPTY = 0x0,
         SD_WAKE = 0x1,
         SD_SLEEP = 0x2,
         SD_REMOVE = 0x4
     };

    StateDialog(JobId jobId, int flags);
    ~StateDialog();

protected:
    enum Selection {
        SD_COMBO_WAKE,
        SD_COMBO_SLEEP,
        SD_COMBO_REMOVE
    };

    void LoadPersistentSettings();
    void SavePersistentSettings();

protected slots:
    void OnOk();

signals:
    void WakeJob(JobId jobId, std::string &password);
    void SleepJob(JobId jobId, std::string &password);
    void RemoveJob(JobId jobId, std::string &password);

private:
    Ui::StateDialog mUi;
    JobId mJobId;
};

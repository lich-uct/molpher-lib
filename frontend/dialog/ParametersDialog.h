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

#include "auxiliary/QtMocHack.h"

#include <QtGui/QDialog>

#include "MolpherParam.h"
#include "global_types.h"
#include "auxiliary/PasswordCache.h"
#include "ui_ParametersDialog.h"

class ParametersDialog : public QDialog
{
    Q_OBJECT

public:
    ParametersDialog(MolpherParam &params, QWidget *parent);
    ParametersDialog(JobId jobId, MolpherParam &params);
    ~ParametersDialog();
    void FillParamsValues(MolpherParam &param);

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();
    void Init(MolpherParam &params);

protected slots:
    void OnButtonDefault();
    void OnOk();

signals:
    void SetParams(MolpherParam &params);
    void SetParams(JobId jobId, MolpherParam &params, std::string &password);

private:
    Ui::ParametersDialog mUi;
    bool mHaveJobId;
    int mJobId;
};

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
#include "ParametersDialog.h"

ParametersDialog::ParametersDialog(MolpherParam &params, QWidget *parent) :
    QDialog(parent),
    mHaveJobId(false)
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);

    Init(params);
}

ParametersDialog::ParametersDialog(JobId jobId, MolpherParam &params) :
    QDialog(GetMainWindow()),
    mHaveJobId(true),
    mJobId(jobId)
{
    mUi.setupUi(this);
    setAttribute(Qt::WA_DeleteOnClose, true);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);

    Init(params);
    connect(this, SIGNAL(SetParams(JobId, MolpherParam &, std::string &)),
        &gCommunicator, SLOT(SetParams(JobId, MolpherParam &, std::string &)));
}

ParametersDialog::~ParametersDialog()
{
    SavePersistentSettings();
}

void ParametersDialog::Init(MolpherParam &params)
{
    FillParamsValues(params);
    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));
    connect(mUi.buttonDefault, SIGNAL(clicked()), this, SLOT(OnButtonDefault()));
    LoadPersistentSettings();
}

void ParametersDialog::FillParamsValues(MolpherParam &params)
{
    mUi.spinCntCandidatesToKeep->setValue(params.cntCandidatesToKeep);
    mUi.spinCntCandidatesToKeepMax->setValue(params.cntCandidatesToKeepMax);
    mUi.spinCntIterations->setValue(params.cntIterations);
    mUi.spinCntMaxMorphs->setValue(params.cntMaxMorphs);
    mUi.spinCntMorphs->setValue(params.cntMorphs);
    mUi.spinCntMorphsInDepth->setValue(params.cntMorphsInDepth);
    mUi.doubleSpinDistToTargetDepthSwitch->setValue(params.distToTargetDepthSwitch);
    mUi.spinItThreshold->setValue(params.itThreshold);
    long int time = params.timeMaxSeconds;
    mUi.spinTimeMaxHours->setValue(time / 3600);
    time = time % 3600;
    mUi.spinTimeMaxMinutes->setValue(time / 60);
    mUi.spinTimeMaxSeconds->setValue(time % 60);
    mUi.doubleSpinMinAcceptableMolecularWeight->setValue(
        params.minAcceptableMolecularWeight);
    mUi.doubleSpinMaxAcceptableMolecularWeight->setValue(
        params.maxAcceptableMolecularWeight);
    
    mUi.chbUseSyntetizedFilter->setChecked(params.useSyntetizedFeasibility);
}

void ParametersDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ParametersDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void ParametersDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("ParametersDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void ParametersDialog::OnButtonDefault()
{
    MolpherParam param;
    FillParamsValues(param);
}

void ParametersDialog::OnOk()
{
    MolpherParam params;
    params.cntCandidatesToKeep = mUi.spinCntCandidatesToKeep->value();
    params.cntCandidatesToKeepMax = mUi.spinCntCandidatesToKeepMax->value();
    params.cntIterations = mUi.spinCntIterations->value();
    params.cntMaxMorphs = mUi.spinCntMaxMorphs->value();
    params.cntMorphs = mUi.spinCntMorphs->value();
    params.cntMorphsInDepth = mUi.spinCntMorphsInDepth->value();
    params.distToTargetDepthSwitch = mUi.doubleSpinDistToTargetDepthSwitch->value();
    params.itThreshold = mUi.spinItThreshold->value();
    params.timeMaxSeconds = mUi.spinTimeMaxSeconds->value() +
        mUi.spinTimeMaxMinutes->value() * 60 + mUi.spinTimeMaxHours->value() * 3600;
    params.minAcceptableMolecularWeight =
        mUi.doubleSpinMinAcceptableMolecularWeight->value();
    params.maxAcceptableMolecularWeight =
        mUi.doubleSpinMaxAcceptableMolecularWeight->value();
    params.useSyntetizedFeasibility = 
            mUi.chbUseSyntetizedFilter->isChecked();
    if (mHaveJobId) {
        std::string password;
        if (PasswordCache::ResolvePassword(mJobId, password)) {
            emit SetParams(mJobId, params, password);
        } else {
            ShowPermissionError();
        }
    } else {
        emit SetParams(params);
    }
}

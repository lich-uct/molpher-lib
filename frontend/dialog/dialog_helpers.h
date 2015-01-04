/*
 Copyright (c) 2012 Martin Straka
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

#pragma once

#include <QtGui/QMainWindow>
#include <QtGui/QFileDialog>
#include <QtGui/QComboBox>
#include <QtGui/QListWidget>

#include "meta/MolpherMoleculeMeta.h"

QString GetFileNameReadFile(QWidget *caller);
QString GetFileNameWriteSDF(QWidget *caller);
QString GetFileNameReadSNP(QWidget *caller);
QString GetFileNameWriteSNP(QWidget *caller);
QStringList GetFileNamesReadSNP(QWidget *caller);

void ShowPermissionError();
void ShowWarning(const QString &text);
void ShowInformation(const QString &text);

void FillListWithSmiles(QListWidget *list, std::vector<MolpherMolecule> &molecules);
void FillListWithFormulas(QListWidget *list, std::vector<MolpherMolecule> &molecules);
void FillListWithChemOperSelectors(std::vector<boost::int32_t> *list);

void FillComboWithFingerprints(QComboBox *comboBox, bool extended = true);
void FillComboWithSimCoeffs(QComboBox *comboBox);
void FillComboWithDimRed(QComboBox *comboBox);

void SelectAll(QListWidget *list);

QMainWindow *GetMainWindow();

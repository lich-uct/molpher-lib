/*
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

#include <string>
#include <vector>
#include <map>
#include <deque>

#include <QtGui/QDialog>

#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"

#include "global_types.h"
#include "MolpherMolecule.h"
#include "ui_FilterDialog.h"

class FilterDialog : public QDialog
{
    Q_OBJECT

public:
    FilterDialog();
    ~FilterDialog();

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();
    void UpdateSelection();

    Fingerprint *GetFingerprint(FingerprintSelector &fingerprintSel, RDKit::ROMol &mol);
    double GetSimCoeff(SimCoeffSelector &simcoeffSel, Fingerprint *fp1, Fingerprint *fp2);

protected slots:
    void OnSearchStructure();
    void OnSearchSimilarity();
    void OnAddSelected();
    void OnAddBookmarks();
    void OnStepBack();
    void OnSaveResults();
    void OnClose();

signals:
    void GetSelectedMolecules(std::vector<MolpherMolecule> &molecules);

private:
    typedef std::map<std::string, MolpherMolecule> FilterSelection;
    typedef std::deque<FilterSelection> FilterHistory;

    Ui::FilterDialog mUi;
    FilterSelection mCurrentSelection;
    FilterHistory mSearchHistory;
};

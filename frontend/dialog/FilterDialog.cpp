/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Peter Szepe

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

#include <cassert>
#include <vector>

#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SanitException.h>
#include <DataStructs/BitOps.h>

#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"

#include "auxiliary/GlobalObjectsHolder.h"
#include "dialog_helpers.h"
#include "inout.h"
#include "FilterDialog.h"

const double simCoeffDefaultValue = 0.8;

FilterDialog::FilterDialog() : QDialog(GetMainWindow())
{
    mUi.setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    setAttribute(Qt::WA_DeleteOnClose, true);

    FillComboWithFingerprints(mUi.comboBoxFingerprint, false);
    FillComboWithSimCoeffs(mUi.comboBoxSimCoeff);
    LoadPersistentSettings();

    connect(mUi.pushButtonStructure, SIGNAL(clicked()), this, SLOT(OnSearchStructure()));
    connect(mUi.pushButtonSimilarity, SIGNAL(clicked()), this, SLOT(OnSearchSimilarity()));
    connect(mUi.pushButtonSelected, SIGNAL(clicked()), this, SLOT(OnAddSelected()));
    connect(mUi.pushButtonBookmarks, SIGNAL(clicked()), this, SLOT(OnAddBookmarks()));
    connect(mUi.pushButtonUndo, SIGNAL(clicked()), this, SLOT(OnStepBack()));
    connect(mUi.pushButtonSave, SIGNAL(clicked()), this, SLOT(OnSaveResults()));
    connect(mUi.pushButtonClose, SIGNAL(clicked()), this, SLOT(OnClose()));
}

FilterDialog::~FilterDialog()
{
    SavePersistentSettings();
}

void FilterDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("FilterDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }

    mUi.comboBoxFingerprint->setCurrentIndex(
        settings.value("fingerprintSel", DEFAULT_FP).toInt());
    mUi.comboBoxSimCoeff->setCurrentIndex(
        settings.value("simCoeffSel", DEFAULT_SC).toInt());
    mUi.doubleSpinBoxSimCoeff->setValue(
        settings.value("simCoeffSelValue", simCoeffDefaultValue ).toDouble());
    mUi.lineEditQuery->setText(
        settings.value("query").toString());
    settings.endGroup();
}

void FilterDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("FilterDialog");
    settings.setValue("geometry", geometry());
    settings.setValue("fingerpringSel", mUi.comboBoxFingerprint->currentIndex());
    settings.setValue("simCoeffSel", mUi.comboBoxSimCoeff->currentIndex());
    settings.setValue("simCoeffSelValue", mUi.doubleSpinBoxSimCoeff->value());
    settings.setValue("query", mUi.lineEditQuery->text());
    settings.endGroup();
}

void FilterDialog::UpdateSelection()
{
    mUi.listWidgetMolecules->clear();
    FilterSelection::iterator it;
    for (it = mCurrentSelection.begin(); it != mCurrentSelection.end(); ++it) {
        mUi.listWidgetMolecules->addItem(QString::fromStdString(it->second.smile));
    }
}

Fingerprint *FilterDialog::GetFingerprint(FingerprintSelector &fingerprintSel, RDKit::ROMol &mol)
{
    Fingerprint *result;
    switch(fingerprintSel) {
    case FP_TOPOLOGICAL:
        result = RDKit::RDKFingerprintMol(mol);
        break;
    case FP_TOPOLOGICAL_LAYERED_1:
        result = RDKit::LayeredFingerprintMol(mol);
        break;
    case FP_TOPOLOGICAL_LAYERED_2:
        result = RDKit::PatternFingerprintMol(mol);
        break;
    case FP_ATOM_PAIRS:
        result = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(mol);
        break;
    case FP_TOPOLOGICAL_TORSION:
        result = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(mol);
        break;
    case FP_MORGAN:
        result = RDKit::MorganFingerprints::getFingerprintAsBitVect(mol, 1, 2048);
        break;
    default:
        assert(false);
        break;
    }

    return result;
}

double FilterDialog::GetSimCoeff(SimCoeffSelector &simcoeffSel, Fingerprint *fp1, Fingerprint *fp2)
{
    double result;

    switch (simcoeffSel) {
    case SC_TANIMOTO:
        result = TanimotoSimilarity(*fp1, *fp2);
        break;
    case SC_COSINE:
        result = CosineSimilarity(*fp1, *fp2);
        break;
    case SC_KULCZYNSKI:
        result = KulczynskiSimilarity(*fp1, *fp2);
        break;
    case SC_DICE:
        result = DiceSimilarity(*fp1, *fp2);
        break;
    case SC_TVERSKY_SUPERSTRUCTURE:
        /* Setting the weighting of prototype (in our case either target or
           neighborhood center) features to 90% (a=0.9) and variant (i.e. morph)
           features to 10% (b=0.1) means that mainly the prototype features
           are important, i.e., this produces a "superstucture-likeness" measure.
           In this case, a Tversky similarity value of 1.0 means that almost all
           prototype features are represented in the variant, 0.0 that almost
           none are. */
        result = TverskySimilarity(*fp1, *fp2, 0.9, 0.1);
        break;
    case SC_TVERSKY_SUBSTRUCTURE:
        /* Conversely, setting the weights to 10% prototype (a=0.1) and
           90% variant (b=0.9) produces a "substructure-likeness" measure,
           where variants almost completely embedded into prototype have
           have values near 1.0. */
        result = TverskySimilarity(*fp1, *fp2, 0.1, 0.9);
        break;
    case SC_SOKAL:
        result = SokalSimilarity(*fp1, *fp2);
        break;
    case SC_MC_CONNAUGHEY:
        result = McConnaugheySimilarity(*fp1, *fp2);
        break;
    case SC_ASYMMETRIC:
        result = AsymmetricSimilarity(*fp1, *fp2);
        break;
    case SC_BRAUN_BLANQUET:
        result = BraunBlanquetSimilarity(*fp1, *fp2);
        break;
    case SC_RUSSEL:
        result = RusselSimilarity(*fp1, *fp2);
        break;
    case SC_ON_BIT:
        result = OnBitSimilarity(*fp1, *fp2);
        break;
    case SC_ALL_BIT:
        result = AllBitSimilarity(*fp1, *fp2);
        break;
    default:
        assert(false);
        break;
    }

    return result;
}

void FilterDialog::OnSearchStructure()
{
    std::string smarts = mUi.lineEditQuery->text().toStdString();

    RDKit::RWMol *patternMol = NULL;
    try {
        patternMol = RDKit::SmartsToMol(smarts);
        if (patternMol) {
            RDKit::MolOps::Kekulize(*patternMol);
        }
    } catch (const ValueErrorException &exc) {
        SynchCout("Cannot kekulize output molecule.");
        delete patternMol;
        patternMol = NULL;
    } catch (const RDKit::MolSanitizeException &exc) {
        SynchCout("MolSanitizeException");
        delete patternMol;
        patternMol = NULL;
    } catch (const std::exception &exc) {
        SynchCout("std::exception");
        delete patternMol;
        patternMol = NULL;
    }

    if (!patternMol) {
        ShowWarning(tr("The SMARTS is invalid."));
        return;
    }

    FilterSelection newSelection;
    FilterSelection::iterator it;

    RDKit::MatchVectType machVect;
    for (it = mCurrentSelection.begin(); it != mCurrentSelection.end(); ++it) {
        RDKit::RWMol *mol = NULL;
        try {
            mol = RDKit::SmilesToMol(it->first);
            if (mol) {
                RDKit::MolOps::Kekulize(*mol);
            } else {
                throw ValueErrorException("");
            }
        } catch (const ValueErrorException &exc) {
            //SynchCout("Cannot kekulize output molecule.");
        }
        if (!mol) {
            continue;
        }

        if (RDKit::SubstructMatch(*mol, *patternMol, machVect)) {
            newSelection.insert(*it);
        }

        delete mol;
    }
    delete patternMol;

    mSearchHistory.push_back(mCurrentSelection);
    mCurrentSelection = newSelection;

    UpdateSelection();
}

void FilterDialog::OnSearchSimilarity()
{
    std::string smiles = mUi.lineEditQuery->text().toStdString();
    FingerprintSelector fingerprintSel =
        (FingerprintSelector) mUi.comboBoxFingerprint->currentIndex();
    SimCoeffSelector simcoeffSel =
        (SimCoeffSelector) mUi.comboBoxSimCoeff->currentIndex();
    double simcoeffLimit = mUi.doubleSpinBoxSimCoeff->value();

    RDKit::RWMol *patternMol = NULL;
    try {
        patternMol = RDKit::SmilesToMol(smiles);
        if (patternMol) {
            RDKit::MolOps::Kekulize(*patternMol);
        }
    } catch (const ValueErrorException &exc) {
        SynchCout("Cannot kekulize output molecule.");
        delete patternMol;
        patternMol = NULL;
    } catch (const RDKit::MolSanitizeException &exc) {
        SynchCout("MolSanitizeException");
        delete patternMol;
        patternMol = NULL;
    } catch (const std::exception &exc) {
        SynchCout("std::exception");
        delete patternMol;
        patternMol = NULL;
    }
    if (!patternMol) {
        ShowWarning(tr("The SMILES is invalid."));
        return;
    }

    Fingerprint *patternFp = GetFingerprint(fingerprintSel, *patternMol);

    FilterSelection newSelection;
    FilterSelection::iterator it;
    for (it = mCurrentSelection.begin(); it != mCurrentSelection.end(); ++it) {
        RDKit::RWMol *mol = NULL;
        try {
            mol = RDKit::SmilesToMol(it->first);
            if (mol) {
                RDKit::MolOps::Kekulize(*mol);
            } else {
                throw ValueErrorException("");
            }
        } catch (const ValueErrorException &exc) {
            //SynchCout("Cannot kekulize output molecule.");
        }
        if (!mol) {
            continue;
        }

        Fingerprint *fp = GetFingerprint(fingerprintSel, *mol);

        double simcoeff = GetSimCoeff(simcoeffSel, patternFp, fp);
        if (simcoeff >= simcoeffLimit) {
            newSelection.insert(*it);
        }

        delete mol;
        delete fp;
    }
    delete patternMol;
    delete patternFp;

    mSearchHistory.push_back(mCurrentSelection);
    mCurrentSelection = newSelection;

    UpdateSelection();
}

void FilterDialog::OnAddSelected()
{
    std::vector<MolpherMolecule> selected;
    emit GetSelectedMolecules(selected);

    std::vector<MolpherMolecule>::iterator it;
    for (it = selected.begin(); it != selected.end(); ++it) {
        mCurrentSelection[it->smile] = *it;
    }

    UpdateSelection();
}

void FilterDialog::OnAddBookmarks()
{
    std::vector<MolpherMolecule> bookmarks;
    gGlobalObjectsHolder.GetBookmarks()->GetCurrent(bookmarks);

    std::vector<MolpherMolecule>::iterator it;
    for (it = bookmarks.begin(); it != bookmarks.end(); ++it) {
        mCurrentSelection[it->smile] = *it;
    }

    UpdateSelection();
}

void FilterDialog::OnStepBack()
{
    if (!mSearchHistory.empty()) {
        mCurrentSelection = mSearchHistory.back();
        mSearchHistory.pop_back();
        UpdateSelection();
    }
}

void FilterDialog::OnSaveResults()
{
    QString filename = GetFileNameWriteSDF(this);
    WriteMolphMolsToSDF(filename.toStdString(), mCurrentSelection);
}

void FilterDialog::OnClose()
{
    this->reject();
}

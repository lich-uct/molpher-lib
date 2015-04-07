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

#include <cassert>

#include <QtCore/QSettings>
#include <QtCore/QObject>
#include <QtGui/QApplication>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>

#include "chemoper_selectors.h"
#include "dialog_helpers.h"
#include "dimred_selectors.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"

QString GetFileNameReadFile(QWidget *caller)
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DialogHelpers");

    QFileDialog::Options options;
    options |= QFileDialog::DontUseNativeDialog;
    QString selectedFilter;
    QFileInfo file(settings.value("loadFile").toString());
    QString fileName = QFileDialog::getOpenFileName(caller,
        QObject::tr("Choose file with molecules"),
        file.dir().path(),
        QObject::tr("SDF Files (*.sdf);;MDL Files (*.mol);;TXT Files (*.txt);;All Files (*)"),
        &selectedFilter,
        options);

    if (!fileName.isEmpty()) {
        settings.setValue("loadFile", fileName);
    }
    settings.endGroup();

    return fileName;
}

QString GetFileNameWriteSDF(QWidget *caller)
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DialogHelpers");

    QFileDialog::Options options;
    options |= QFileDialog::DontUseNativeDialog;
    QString selectedFilter;
    QFileInfo file(settings.value("saveSDF").toString());
    QString fileName = QFileDialog::getSaveFileName(caller,
        QObject::tr("Save as SDF file"),
        file.dir().path(),
        QObject::tr("SDF Files (*.sdf);;All Files (*)"),
        &selectedFilter,
        options);

    if (!fileName.isEmpty()) {
        QString suffix(".sdf");
        QRegExp reg("\\.[sS][dD][fF]$");
        fileName.replace(reg, suffix);

        int index = fileName.lastIndexOf(suffix);
        if (index != fileName.length() - suffix.length()) {
            fileName += suffix;
        }
        settings.setValue("saveSDF", fileName);
    }
    settings.endGroup();

    return fileName;
}

QString GetFileNameReadSNP(QWidget *caller)
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DialogHelpers");

    QFileDialog::Options options;
    options |= QFileDialog::DontUseNativeDialog;
    QString selectedFilter;
    QFileInfo file(settings.value("loadSNP").toString());
    QString fileName = QFileDialog::getOpenFileName(caller,
        QObject::tr("Choose SNP file"),
        file.dir().path(),
        QObject::tr("SNP Files (*.snp);;All Files (*)"),
        &selectedFilter,
        options);

    if (!fileName.isEmpty()) {
        settings.setValue("loadSNP", fileName);
    }
    settings.endGroup();

    return fileName;
}

QString GetFileNameWriteSNP(QWidget *caller)
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DialogHelpers");

    QFileDialog::Options options;
    options |= QFileDialog::DontUseNativeDialog;
    QString selectedFilter;
    QFileInfo file(settings.value("saveSNP").toString());
    QString fileName = QFileDialog::getSaveFileName(caller,
        QObject::tr("Save as SNP file"),
        file.dir().path(),
        QObject::tr("SNP Files (*.snp);;All Files (*)"),
        &selectedFilter,
        options);

    if (!fileName.isEmpty()) {
        QString suffix(".snp");
        QRegExp reg("\\.[sS][nN][pP]$");
        fileName.replace(reg, suffix);

        int index = fileName.lastIndexOf(suffix);
        if (index != fileName.length() - suffix.length()) {
            fileName += suffix;
        }
        settings.setValue("saveSNP", fileName);
    }
    settings.endGroup();

    return fileName;
}

QStringList GetFileNamesReadSNP(QWidget *caller)
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DialogHelpers");

    QFileDialog::Options options;
    options |= QFileDialog::DontUseNativeDialog;
    QString selectedFilter;
    QStringList filenames = QFileDialog::getOpenFileNames(caller,
        QObject::tr("Choose SNP files"),
        QObject::tr(settings.value("loadSNP").toByteArray()),
        QObject::tr("SNP Files (*.snp);;All Files (*)"),
        &selectedFilter,
        options);

    if (!filenames.isEmpty() && !filenames.first().isEmpty()) {
        QFileInfo file(filenames.first());
        settings.setValue("loadSNP", file.dir().path());
    }
    settings.endGroup();

    return filenames;
}

void ShowPermissionError()
{
    QMessageBox msgBox(GetMainWindow());
    msgBox.setText(QObject::tr("You don't have write permissions for this job."));
    msgBox.setStandardButtons(QMessageBox::Discard);
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.exec();
}

void ShowWarning(const QString &text)
{
    QMessageBox msgBox(GetMainWindow());
    msgBox.setText(text);
    msgBox.setStandardButtons(QMessageBox::Discard);
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.exec();
}

void ShowInformation(const QString &text)
{
    QMessageBox msgBox(GetMainWindow());
    msgBox.setText(text);
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();
}

void FillListWithSmiles(QListWidget *list, std::vector<MolpherMolecule> &molecules)
{
    QString qSmiles;
    std::vector<MolpherMolecule>::iterator it;
    for (it = molecules.begin(); it != molecules.end(); ++it) {
        qSmiles = QString(it->smile.c_str());
        QListWidgetItem *item = new QListWidgetItem(qSmiles, list);
        list->addItem(item);
    }
}

void FillListWithFormulas(QListWidget *list, std::vector<MolpherMolecule> &molecules)
{
    QString qFormula;
    std::vector<MolpherMolecule>::iterator it;
    for (it = molecules.begin(); it != molecules.end(); ++it) {
        qFormula = QString(it->formula.c_str());
        QListWidgetItem *item = new QListWidgetItem(qFormula, list);
        list->addItem(item);
    }
}

void FillListWithChemOperSelectors(std::vector<boost::int32_t> *list)
{
    list->push_back(OP_ADD_ATOM);
    list->push_back(OP_REMOVE_ATOM);
    list->push_back(OP_ADD_BOND);
    list->push_back(OP_REMOVE_BOND);
    list->push_back(OP_MUTATE_ATOM);
    list->push_back(OP_INTERLAY_ATOM);
    list->push_back(OP_BOND_REROUTE);
    list->push_back(OP_BOND_CONTRACTION);
}

void FillComboWithFingerprints(QComboBox *comboBox, bool extended)
{
    QStringList listFingerprints;
    listFingerprints.append(FingerprintLongDesc(FP_ATOM_PAIRS));
    listFingerprints.append(FingerprintLongDesc(FP_MORGAN));
    listFingerprints.append(FingerprintLongDesc(FP_TOPOLOGICAL));
    listFingerprints.append(FingerprintLongDesc(FP_TOPOLOGICAL_LAYERED_1));
    listFingerprints.append(FingerprintLongDesc(FP_TOPOLOGICAL_LAYERED_2));
    listFingerprints.append(FingerprintLongDesc(FP_TOPOLOGICAL_TORSION));
    if (extended) {
        listFingerprints.append(FingerprintLongDesc(FP_EXT_ATOM_PAIRS));
        listFingerprints.append(FingerprintLongDesc(FP_EXT_MORGAN));
        listFingerprints.append(FingerprintLongDesc(FP_EXT_TOPOLOGICAL));
        listFingerprints.append(FingerprintLongDesc(FP_EXT_TOPOLOGICAL_LAYERED_1));
        listFingerprints.append(FingerprintLongDesc(FP_EXT_TOPOLOGICAL_LAYERED_2));
        listFingerprints.append(FingerprintLongDesc(FP_EXT_TOPOLOGICAL_TORSION));
    }
    listFingerprints.append(FingerprintLongDesc(FP_VECTORFP));
    comboBox->addItems(listFingerprints);
    comboBox->setCurrentIndex(DEFAULT_FP);
}

void FillComboWithSimCoeffs(QComboBox *comboBox)
{
    QStringList listSimCoeffs;
    listSimCoeffs.append(SimCoeffLongDesc(SC_ALL_BIT));
    listSimCoeffs.append(SimCoeffLongDesc(SC_ASYMMETRIC));
    listSimCoeffs.append(SimCoeffLongDesc(SC_BRAUN_BLANQUET));
    listSimCoeffs.append(SimCoeffLongDesc(SC_COSINE));
    listSimCoeffs.append(SimCoeffLongDesc(SC_DICE));
    listSimCoeffs.append(SimCoeffLongDesc(SC_KULCZYNSKI));
    listSimCoeffs.append(SimCoeffLongDesc(SC_MC_CONNAUGHEY));
    listSimCoeffs.append(SimCoeffLongDesc(SC_ON_BIT));
    listSimCoeffs.append(SimCoeffLongDesc(SC_RUSSEL));
    listSimCoeffs.append(SimCoeffLongDesc(SC_SOKAL));
    listSimCoeffs.append(SimCoeffLongDesc(SC_TANIMOTO));
    listSimCoeffs.append(SimCoeffLongDesc(SC_TVERSKY_SUBSTRUCTURE));
    listSimCoeffs.append(SimCoeffLongDesc(SC_TVERSKY_SUPERSTRUCTURE));
    comboBox->addItems(listSimCoeffs);
    comboBox->setCurrentIndex(DEFAULT_SC);
}

void FillComboWithDimRed(QComboBox *comboBox)
{
    QStringList listDimRed;
    listDimRed.append(DimRedLongDesc(DR_KAMADAKAWAI));
    listDimRed.append(DimRedLongDesc(DR_PCA));    
    comboBox->addItems(listDimRed);
    comboBox->setCurrentIndex(DEFAULT_DR);
}

void SelectAll(QListWidget *list)
{
    list->selectAll();
    list->setFocus();
}

QMainWindow *GetMainWindow()
{
    QMainWindow *main = NULL;
    foreach(QWidget *widget, QApplication::topLevelWidgets()) {
        if(widget->objectName() == "MainWindow") {
            main = qobject_cast<QMainWindow *>(widget);
        }
    }
    assert(main);
    return main;
}

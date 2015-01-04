/*
 Copyright (c) 2012 Marek Mikes
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

#include <boost/date_time/posix_time/posix_time.hpp>

#include <cassert>

#include <QtCore/QString>
#include <QtCore/QByteArray>
#include <QtCore/QBuffer>
#include <QtCore/QUrl>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QSettings>
#include <QtCore/QProcess>
#include <QtGui/QToolBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QPushButton>
#include <QtGui/QInputDialog>
#include <QtGui/QIcon>
#include <QtGui/QImage>
#include <QtGui/QPixmap>
#include <QtGui/QApplication>

#include "FrontendCommunicator.h"
#include "meta/MolpherMoleculeMeta.h"
#include "auxiliary/GlobalObjectsHolder.h"
#include "Bookmarks.h"
#include "dialog/dialog_helpers.h"
#include "MoleculePainter.h"

/**
 * Return image for given molecule. Use extern program indigo-depict to render
 * the image molecule's structure. The image is for short time store in 
 * frontend's Temp folder.
 * @param molecule
 * @return 
 */
QImage renderMolecule(const MolpherMolecule &molecule)
{
    std::string filename;
    filename += QDir::currentPath().toStdString();
    filename += "/Temp/";
    filename += gFrontendId + "/";
    filename += boost::posix_time::to_iso_string(
        boost::posix_time::second_clock::local_time());
    filename += ".png";    
    
    QString program("indigo-depict");
    
    QStringList arguments;
    arguments << "-";
    arguments << QString::fromStdString(molecule.smile);
    arguments << QString::fromStdString(filename);
    arguments << "-bond" << "25";
    arguments << "-comment" << QString::fromStdString(molecule.formula);
    arguments << "-commentoffset" << "8";
    arguments << "-commentsize" << "14";
    arguments << "-commentalign" << "0";
    
    QImage image;
    if (QProcess::execute(program, arguments) == 0) {
        image = QImage(QString::fromStdString(filename));
        QFile::remove(QString::fromStdString(filename));
    } else {
        image = MoleculePainter::DrawMolecule(molecule.smile).toImage();
    }
    
    return image;
}

/**
 * Add image to the given item based on given molecule representation.
 * @param item
 * @param molecule
 */
void updateItem(QListWidgetItem& item, const MolpherMolecule &molecule)
{
    QImage image = renderMolecule(molecule);
    
    QPixmap pixmap = QPixmap::fromImage(image);
    QIcon icon;
    item.setIcon(QIcon(pixmap));
}

Bookmarks::Bookmarks(QWidget *parent) :
    QDockWidget(tr("Bookmarks"), parent),
    mDisplayFormula(true)
{
    setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    mToolBoxBookmarks = new QToolBox(this);
    LoadPersistentSettings();

    QPushButton *buttonAddGroup = new QPushButton(tr("New group"), this);
    connect(buttonAddGroup, SIGNAL(clicked()), this, SLOT(OnButtonAddGroup()));

    QPushButton *buttonDeleteGroup = new QPushButton(tr("Delete group"), this);
    connect(buttonDeleteGroup, SIGNAL(clicked()), this, SLOT(OnButtonDeleteGroup()));

    QVBoxLayout *mainLayout = new QVBoxLayout();
    QHBoxLayout *buttonLayout = new QHBoxLayout();
    mainLayout->addWidget(mToolBoxBookmarks);
    buttonLayout->addWidget(buttonAddGroup);
    buttonLayout->addWidget(buttonDeleteGroup);
    mainLayout->addLayout(buttonLayout);
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);

    QWidget *widget = new QWidget(this);
    widget->setLayout(mainLayout);
    this->setWidget(widget);

    gGlobalObjectsHolder.Add(this);
}

Bookmarks::~Bookmarks()
{
    gGlobalObjectsHolder.Remove(this);
    SavePersistentSettings();
}

void Bookmarks::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Alt:
        mDisplayFormula = !mDisplayFormula;
        UpdateMolLabels();
        break;
    case Qt::Key_Delete:
        DeleteSelectedBookmarks();
        break;
    default:
        break;
    }
    QDockWidget::keyPressEvent(event);
}

void Bookmarks::closeEvent(QCloseEvent *event)
{
    emit CloseBookmarks();
}

void Bookmarks::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("Bookmarks");
    QStringList listGroup(settings.childGroups());
    for (int i = 0; i < listGroup.count(); ++i) {
        QListWidget *newGroup = new QListWidget(this);
        newGroup->setSelectionMode(QAbstractItemView::ExtendedSelection);
        mToolBoxBookmarks->addItem(newGroup, listGroup[i]);

        settings.beginGroup(listGroup[i]);
        QStringList listKeys(settings.childKeys());
        for (int j = 0; j < listKeys.count(); ++j) {
            std::string formula = settings.value(listKeys[j]).toString().toStdString();
            std::string smile = listKeys[j].toStdString();
            if (formula.size() > 0) {
                MolpherMolecule mol(smile, formula);

                QListWidget *choosenGroup = qobject_cast<QListWidget *>(
                    mToolBoxBookmarks->widget(i));
                choosenGroup->setIconSize(QSize(100,100));
                
                QString label;
                if (mDisplayFormula) {
                    label = QString::fromStdString(formula);
                } else {
                    label = QString::fromStdString(smile);
                }
                QVariant variant;
                QListWidgetItem *item = new QListWidgetItem(label, choosenGroup);
                variant = QVariant::fromValue(mol);
                item->setData(Qt::UserRole, variant);
                updateItem(*item, mol);
                choosenGroup->addItem(item);
            }
        }
        settings.endGroup();
    }
    if (mToolBoxBookmarks->count() > 0) {
        mToolBoxBookmarks->setCurrentIndex(settings.value("currentIndex", 0).toInt());
    }
    settings.endGroup();
}

void Bookmarks::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("Bookmarks");
    settings.setValue("currentIndex", mToolBoxBookmarks->currentIndex());
    settings.endGroup();
}



bool Bookmarks::Add(MolpherMolecule &molecule)
{
    int currentIndex = mToolBoxBookmarks->currentIndex();
    if (currentIndex >= 0) {
        QListWidget *currentGroup = qobject_cast<QListWidget *>(
            mToolBoxBookmarks->widget(currentIndex));

        QString label;
        if (mDisplayFormula) {
            label = QString::fromStdString(molecule.formula);
        } else {
            label = QString::fromStdString(molecule.smile);
        }
        QVariant variant;
        QListWidgetItem *item = new QListWidgetItem(label, currentGroup);
        variant = QVariant::fromValue(molecule);
        item->setData(Qt::UserRole, variant);
        updateItem(*item, molecule);
        currentGroup->addItem(item);
        
        QSettings settings("siret", "molpher");
        settings.beginGroup("Bookmarks");

        settings.beginGroup(mToolBoxBookmarks->itemText(currentIndex));
        settings.setValue(QString::fromStdString(molecule.smile),
            QString::fromStdString(molecule.formula));
        settings.endGroup();

        settings.endGroup();
        return true;
    } else {
        return false;
    }
}

QStringList Bookmarks::GetAllGroups()
{
    QStringList listGroup;
    for (int i = 0; i < mToolBoxBookmarks->count(); ++i) {
        listGroup.append(mToolBoxBookmarks->itemText(i));
    }
    return listGroup;
}

void Bookmarks::GetAll(std::vector<MolpherMolecule> &mols) const
{
    for (int i = 0; i < mToolBoxBookmarks->count(); ++i) {
        QListWidget *oneGroup = qobject_cast<QListWidget *>(
            mToolBoxBookmarks->widget(i));
        for (int j = 0; j < oneGroup->count(); ++j) {
            mols.push_back(
                oneGroup->item(j)->data(Qt::UserRole).value<MolpherMolecule>());
        }
    }
}

void Bookmarks::GetCurrent(std::vector<MolpherMolecule> &mols) const
{
    QListWidget *currentGroup = qobject_cast<QListWidget *>(
        mToolBoxBookmarks->currentWidget());
    if (currentGroup) {
        for (int j = 0; j < currentGroup->count(); ++j) {
            mols.push_back(
                currentGroup->item(j)->data(Qt::UserRole).value<MolpherMolecule>());
        }
    }
}

void Bookmarks::UpdateMolLabels()
{
    for (int i = 0; i < mToolBoxBookmarks->count(); ++i) {
        QListWidget *oneGroup = qobject_cast<QListWidget *>(
            mToolBoxBookmarks->widget(i));
        for (int j = 0; j < oneGroup->count(); ++j) {
            QString label;
            if (mDisplayFormula) {
                label = QString::fromStdString(
                    oneGroup->item(j)->data(Qt::UserRole).value<MolpherMolecule>().formula);
            } else {
                label = QString::fromStdString(
                    oneGroup->item(j)->data(Qt::UserRole).value<MolpherMolecule>().smile);
            }
            oneGroup->item(j)->setText(label);
        }
    }
}

void Bookmarks::OnButtonAddGroup()
{
    bool ok;
    QString text = QInputDialog::getText(this, tr("New group"),
        tr("Fill in a name of new group:"), QLineEdit::Normal, "", &ok);
    if (ok && !text.isEmpty()) {
        QStringList listGroup(GetAllGroups());
        if (!listGroup.contains(text)) {
            QListWidget *newGroup = new QListWidget(this);
            newGroup->setIconSize(QSize(100,100));
            newGroup->setSelectionMode(QAbstractItemView::ExtendedSelection);
            mToolBoxBookmarks->insertItem(InsertIndex(text), newGroup, text);
            mToolBoxBookmarks->setCurrentIndex(mToolBoxBookmarks->indexOf(newGroup));
                        
            QSettings settings("siret", "molpher");
            settings.beginGroup("Bookmarks");

            settings.beginGroup(text);
            settings.setValue("", "");
            settings.endGroup();

            settings.endGroup();
        } else {
            ShowWarning(tr("Group with the same name already exists."));
        }
    }
}

void Bookmarks::DeleteSelectedBookmarks()
{
    int index = mToolBoxBookmarks->currentIndex();
    if (index >= 0) {
        QListWidget *currentGroup = qobject_cast<QListWidget *>(
            mToolBoxBookmarks->widget(index));
        QList<QListWidgetItem *> selectedItems = currentGroup->selectedItems();

        foreach (QListWidgetItem *item, selectedItems) {
            std::string smile = item->data(Qt::UserRole).value<MolpherMolecule>().smile;
            delete currentGroup->takeItem(currentGroup->row(item));

            QSettings settings("siret", "molpher");
            settings.beginGroup("Bookmarks");

            settings.beginGroup(mToolBoxBookmarks->itemText(index));
            settings.remove(QString::fromStdString(smile));
            settings.endGroup();

            settings.endGroup();
        }
    }
}

void Bookmarks::OnButtonDeleteGroup()
{
    int index = mToolBoxBookmarks->currentIndex();
    if (index >= 0) {
        QListWidget *currentGroup = qobject_cast<QListWidget *>(
            mToolBoxBookmarks->widget(index));

        QSettings settings("siret", "molpher");
        settings.beginGroup("Bookmarks");
        settings.remove(mToolBoxBookmarks->itemText(index));
        settings.endGroup();

        mToolBoxBookmarks->removeItem(index);
        delete currentGroup;
    }
}

int Bookmarks::InsertIndex(QString &text)
{
   for (int i = 0; i < mToolBoxBookmarks->count(); ++i) {
       if (text < mToolBoxBookmarks->itemText(i)) {
           return i;
       }
   }
   return mToolBoxBookmarks->count();
}

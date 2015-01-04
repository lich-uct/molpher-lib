/*
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

#include <vector>

#include <QtCore/QVariant>
#include <QtCore/QSettings>
#include <QtGui/QDesktopWidget>
#include <QtCore/QDebug>

#include "inout.h"
#include "components/Bookmarks.h"
#include "dialog_helpers.h"
#include "auxiliary/GlobalObjectsHolder.h"
#include "auxiliary/PasswordCache.h"
#include "FrontendCommunicator.h"
#include "DecoysDialog.h"

DecoysDialog::DecoysDialog(std::vector<MolpherMolecule> &decoys, QWidget *parent) :
    QDialog(parent),
    mHaveJobId(false),
    mDisplayFormula(true)
{
    BasicAdjustment();
    AppendList(mUi.listDecoys, decoys);
}

DecoysDialog::DecoysDialog(JobId jobId, std::vector<MolpherMolecule> &decoys) :
    QDialog(GetMainWindow()),
    mJobId(jobId),
    mHaveJobId(true),
    mDisplayFormula(true)
{
    setAttribute(Qt::WA_DeleteOnClose, true);

    BasicAdjustment();
    AppendList(mUi.listDecoys, decoys);

    connect(this,
        SIGNAL(SetDecoys(JobId, std::vector<MolpherMolecule> &, std::string &)),
        &gCommunicator,
        SLOT(SetDecoys(JobId, std::vector<MolpherMolecule> &, std::string &)));
}

DecoysDialog::~DecoysDialog()
{
    SavePersistentSettings();
}

void DecoysDialog::keyPressEvent(QKeyEvent *e)
{
    switch (e->key()) {
    case Qt::Key_Alt:
        mDisplayFormula = !mDisplayFormula;
        UpdateList(mUi.listBookmarks);
        UpdateList(mUi.listDecoys);
        UpdateList(mUi.listOpened);
        break;
    case Qt::Key_Delete:
        OnButtonDeleteDecoy();
        break;
    default:
        break;
    }
    QDialog::keyPressEvent(e);
}

void DecoysDialog::BasicAdjustment()
{
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    mUi.setupUi(this);

    LoadBookmarks();
    LoadPersistentSettings();

    connect(mUi.listOpened, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
        this, SLOT(OnListOpenedDoubleClicked(QListWidgetItem *)));
    connect(mUi.listBookmarks, SIGNAL(itemDoubleClicked(QListWidgetItem *)),
        this, SLOT(OnListBookmarksDoubleClicked(QListWidgetItem *)));
    connect(mUi.buttonFromOpened, SIGNAL(pressed()), this, SLOT(OnButtonFromOpened()));
    connect(mUi.buttonFromBookmarks, SIGNAL(pressed()), this, SLOT(OnButtonFromBookmarks()));
    connect(mUi.buttonLoad, SIGNAL(pressed()), this, SLOT(OnButtonLoad()));
    connect(mUi.buttonSelectAllOpened, SIGNAL(pressed()),
        this, SLOT(OnButtonSelectAllOpened()));
    connect(mUi.buttonDeleteDecoy, SIGNAL(pressed()), this, SLOT(OnButtonDeleteDecoy()));
    connect(mUi.buttonSelectAllDecoys, SIGNAL(pressed()),
        this, SLOT(OnButtonSelectAllDecoys()));
    connect(mUi.buttonSelectAllBookmarks, SIGNAL(pressed()),
        this, SLOT(OnButtonSelectAllBookmarks()));
    connect(mUi.buttonLoadBookmarks, SIGNAL(clicked()),
        this, SLOT(OnButtonLoadBookmarks()));
    connect(this, SIGNAL(accepted()), this, SLOT(OnOk()));
}

void DecoysDialog::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DecoysDialog");
    QRect rect = settings.value("geometry").toRect();
    if (rect.isEmpty()) {
        move(QApplication::desktop()->screen()->rect().center() - this->rect().center());
    } else {
        setGeometry(rect);
    }
    settings.endGroup();
}

void DecoysDialog::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("DecoysDialog");
    settings.setValue("geometry", geometry());
    settings.endGroup();
}

void DecoysDialog::LoadBookmarks()
{
    std::vector<MolpherMolecule> bookmarks;
    gGlobalObjectsHolder.GetBookmarks()->GetCurrent(bookmarks);
    AppendList(mUi.listBookmarks, bookmarks);
}

void DecoysDialog::OnButtonLoadBookmarks()
{
    LoadBookmarks();
}

void DecoysDialog::FillDecoysList(std::vector<MolpherMolecule> &decoys)
{
    mUi.listDecoys->clear();
    AppendList(mUi.listDecoys, decoys);
}

void DecoysDialog::AppendList(QListWidget *list, std::vector<MolpherMolecule> &molecules)
{
    QString label;
    QVariant variant;
    std::vector<MolpherMolecule>::iterator it;
    for (it = molecules.begin(); it != molecules.end(); ++it) {
        if (mDisplayFormula) {
            label = QString::fromStdString(it->formula);
        } else {
            label = QString::fromStdString(it->smile);
        }
        QListWidgetItem *item = new QListWidgetItem(label, list);
        variant = QVariant::fromValue(*it);
        item->setData(Qt::UserRole, variant);
        list->addItem(item);
    }
}

void DecoysDialog::UpdateList(QListWidget *list)
{
    for (int i = 0; i < list->count(); ++i) {
        QString label;
        if (mDisplayFormula) {
            label = QString::fromStdString(
                list->item(i)->data(Qt::UserRole).value<MolpherMolecule>().formula);
        } else {
            label = QString::fromStdString(
                list->item(i)->data(Qt::UserRole).value<MolpherMolecule>().smile);
        }
        list->item(i)->setText(label);
    }
}

void DecoysDialog::Move(QListWidget *source, QListWidget *destination)
{
    assert(NULL != source);
    assert(NULL != destination);
    assert(source != destination);

    foreach (QListWidgetItem *item, source->selectedItems()) {
        if (IsIn(item, destination)) {
            // item exists in destination
            continue;
        }
        QListWidgetItem * itemDest =  new QListWidgetItem(*item);
        destination->addItem(itemDest);
    }
}

bool DecoysDialog::IsIn(const QListWidgetItem *item, const QListWidget *list) const
{
    // smiles are compared
    std::string inputSmile;
    std::string listSmile;
    QVariant variant;
    for (int i = 0; i < list->count(); ++i) {
        variant = item->data(Qt::UserRole);
        inputSmile = variant.value<MolpherMolecule>().smile;
        variant = list->item(i)->data(Qt::UserRole);
        listSmile = variant.value<MolpherMolecule>().smile;

        if (inputSmile.compare(listSmile) == 0) {
            return true;
        }
    }

    // not found
    return false;
}

void DecoysDialog::GetDecoysFromList(std::vector<MolpherMolecule> &decoys) const
{
    QListWidgetItem *item;
    QVariant variant;
    for (int i = 0; i < mUi.listDecoys->count(); ++i) {
        item = mUi.listDecoys->item(i);
        variant = item->data(Qt::UserRole);
        decoys.push_back(variant.value<MolpherMolecule>());
    }
}

void DecoysDialog::OnListOpenedDoubleClicked(QListWidgetItem *item)
{
    Q_UNUSED(item)
    Move(mUi.listOpened, mUi.listDecoys);
}

void DecoysDialog::OnListBookmarksDoubleClicked(QListWidgetItem *item)
{
    Q_UNUSED(item)
    Move(mUi.listBookmarks, mUi.listDecoys);
}

void DecoysDialog::OnButtonFromOpened()
{
    Move(mUi.listOpened, mUi.listDecoys);
}

void DecoysDialog::OnButtonFromBookmarks()
{
    Move(mUi.listBookmarks, mUi.listDecoys);
}

void DecoysDialog::OnButtonLoad()
{
    std::vector<MolpherMolecule> mols;
    QString fileName = GetFileNameReadFile(this);
    ReadMolphMolsFromFile(fileName.toStdString(), mols);

    AppendList(mUi.listOpened, mols);
}

void DecoysDialog::OnButtonDeleteDecoy()
{
    foreach (QListWidgetItem *item, mUi.listDecoys->selectedItems()) {
        delete mUi.listDecoys->takeItem(mUi.listDecoys->row(item));
    }
}

void DecoysDialog::OnButtonSelectAllOpened()
{
    SelectAll(mUi.listOpened);
}

void DecoysDialog::OnButtonSelectAllDecoys()
{
    SelectAll(mUi.listDecoys);
}

void DecoysDialog::OnButtonSelectAllBookmarks()
{
    SelectAll(mUi.listBookmarks);
}

void DecoysDialog::OnOk()
{
    std::vector<MolpherMolecule> decoys;
    GetDecoysFromList(decoys);

    if (mHaveJobId) {
        std::string password;
        if (PasswordCache::ResolvePassword(mJobId, password)) {
            emit SetDecoys(mJobId, decoys, password);
        } else {
            ShowPermissionError();
        }
    } else {
        emit SetDecoys(decoys);
    }
}

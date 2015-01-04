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

#pragma once

#include <QtGui/QDockWidget>
#include <QtGui/QListWidget>
#include <QtGui/QKeyEvent>

#include "MolpherMolecule.h"

class QToolBox;
class QStringList;

class Bookmarks : public QDockWidget
{
    Q_OBJECT

public:
    Bookmarks(QWidget *parent);
    ~Bookmarks();

    bool Add(MolpherMolecule &molecule);
    void GetAll(std::vector<MolpherMolecule> &mols) const;
    void GetCurrent(std::vector<MolpherMolecule> &mols) const;

protected:
    void keyPressEvent(QKeyEvent *event);
    void closeEvent(QCloseEvent *event);
    void LoadPersistentSettings();
    void SavePersistentSettings();
    void UpdateMolLabels();
    void DeleteSelectedBookmarks();
    QStringList GetAllGroups();
    int InsertIndex(QString &text);

protected slots:
    void OnButtonAddGroup();
    void OnButtonDeleteGroup();

signals:
    void CloseBookmarks();

private:
    bool mDisplayFormula;
    QToolBox *mToolBoxBookmarks;
};

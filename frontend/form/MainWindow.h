/*
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

#include <QtGui/QMainWindow>
#include <QtGui/QLabel>

#include "global_types.h"
#include "components/JobsQueue.h"
#include "components/Bookmarks.h"

namespace Ui {
class MainWindow;
}

class Tab;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void OnCloseTab(int index);
    void OnOpenLiveTab(JobId jobId, IterationSnapshot &snapshot);
    void OnOpenDetachedTab(JobId jobId, IterationSnapshot &snapshot);
    void OnOpenAdhocTab(IterationSnapshotProxy &snapshot);
    void OnOpenAdhocTab(std::vector<IterationSnapshotProxy> &snapshots);
    void OnJobsUpdated(const JobGroup &jobs);
    void OnUpdateJobHistory(JobId jobId, IterIdx maxIterIdx);
    void OnDisplayOnlineState(const QString &ip);
    void OnDisplayOfflineState();

protected:
    void LoadPersistentSettings();
    void SavePersistentSettings();
    int GetLiveTabIndex(JobId jobId);

protected slots:
    void OnOpenSnapshots();
    void OnCreateJob();
    void OnShowJobsQueue(bool show);
    void OnShowBookmarks(bool show);
    void OnConnect();
    void OnDisconnect();
    void OnAbout();
    void OnLegend();
    void OnCloseBookmarks();
    void OnCloseJobsQueue();
    void OnClearSettings();
signals:
    void DisconnectFromBackend();

private:
    Ui::MainWindow *mUi;
    JobsQueue *mJobsQueue;
    Bookmarks *mBookmarks;
    QLabel mConnectionState;
};

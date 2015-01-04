/*
 Copyright (c) 2012 Martin Straka
 Copyright (c) 2012 Marek Mikes
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

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <QtGui/QLabel>
#include <QtCore/QSettings>
#include <QtCore/QStringList>
#include <QtGui/QAction>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QMessageBox>
#include <QtGui/QPixmap>
#include <QtGui/QIcon>

#include "meta/NeighborhoodTaskMeta.h"
#include "meta/IterationSnapshotMeta.h"

#include "auxiliary/GlobalObjectsHolder.h"
#include "components/VisualizedMolecule.h"

#include "inout.h"
#include "dialog/dialog_helpers.h"
#include "FrontendCommunicator.h"
#include "components/JobsQueue.h"
#include "components/Bookmarks.h"
#include "components/Tab.h"
#include "components/ChemicalSpaceView.h"

#include "dialog/CreateSnapshotDialog.h"
#include "dialog/DecoysDialog.h"
#include "dialog/ConnectDialog.h"

#include "ui_MainWindow.h"
#include "MainWindow.h"

#include "Version.hpp"
#include "revision.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    mUi(new Ui::MainWindow)
{
    mUi->setupUi(this);
    setWindowIcon(QIcon(":/icons/molpher.png"));

    mConnectionState.setText(tr("Offline"));
    mUi->actionConnect->setEnabled(true);
    mUi->actionDisconnect->setEnabled(false);
    mUi->actionNewJob->setEnabled(false);
    mUi->statusBar->addPermanentWidget(&mConnectionState);

    mJobsQueue = new JobsQueue(this);
    mJobsQueue->setObjectName("jobsQueue");
    addDockWidget(Qt::BottomDockWidgetArea, mJobsQueue);
    mJobsQueue->hide();
    qRegisterMetaType<MolpherMolecule>();

    mUi->tabWidget->setTabsClosable(true);
    mUi->tabWidget->setMovable(true);

    mBookmarks = new Bookmarks(this);
    mBookmarks->setObjectName("bookmarks");
    addDockWidget(Qt::RightDockWidgetArea, mBookmarks);
    mBookmarks->hide();

    connect(mBookmarks, SIGNAL(CloseBookmarks()), this, SLOT(OnCloseBookmarks()));
    connect(mJobsQueue, SIGNAL(CloseJobsQueue()), this, SLOT(OnCloseJobsQueue()));

    connect(mJobsQueue, SIGNAL(OpenLiveTab(JobId, IterationSnapshot &)),
        this, SLOT(OnOpenLiveTab(JobId, IterationSnapshot &)));
    connect(mJobsQueue, SIGNAL(OpenDetachedTab(JobId, IterationSnapshot &)),
        this, SLOT(OnOpenDetachedTab(JobId, IterationSnapshot &)));
    connect(mJobsQueue, SIGNAL(JobsUpdated(const JobGroup &)),
        this, SLOT(OnJobsUpdated(const JobGroup &)));
    connect(mJobsQueue, SIGNAL(UpdateJobHistory(JobId, IterIdx)),
        this, SLOT(OnUpdateJobHistory(JobId, IterIdx)));

    connect(mUi->tabWidget, SIGNAL(tabCloseRequested(int)), this, SLOT(OnCloseTab(int)));
    connect(mUi->actionOpenSnapshots, SIGNAL(triggered()), this,  SLOT(OnOpenSnapshots()));
    connect(mUi->actionNewJob, SIGNAL(triggered()), this,  SLOT(OnCreateJob()));
    connect(mUi->actionEditJobQueue, SIGNAL(toggled(bool)), this, SLOT(OnShowJobsQueue(bool)));
    connect(mUi->actionBookmarks, SIGNAL(toggled(bool)), this, SLOT(OnShowBookmarks(bool)));

    connect(mUi->actionConnect, SIGNAL(triggered()), this, SLOT(OnConnect()));
    connect(mUi->actionDisconnect, SIGNAL(triggered()), this, SLOT(OnDisconnect()));
    connect(this, SIGNAL(DisconnectFromBackend()),
        &gCommunicator, SLOT(DisconnectFromBackend()));

    connect(&gCommunicator, SIGNAL(DisplayOnlineState(const QString &)),
        this, SLOT(OnDisplayOnlineState(const QString &)), Qt::QueuedConnection);
    connect(&gCommunicator, SIGNAL(DisplayOfflineState()),
        this, SLOT(OnDisplayOfflineState()), Qt::QueuedConnection);

    connect(mUi->actionAbout, SIGNAL(triggered()), this, SLOT(OnAbout()));
    connect(mUi->actionLegend, SIGNAL(triggered()), this, SLOT(OnLegend()));
    connect(mUi->actionClearSettings, SIGNAL(triggered()), this, SLOT(OnClearSettings()));

    qRegisterMetaType<IterationSnapshot>();
    qRegisterMetaType<IterationSnapshotProxy>();

    LoadPersistentSettings();

    if (!mBookmarks->isHidden()) {
        mUi->actionBookmarks->setChecked(true);
    }
    if (!mJobsQueue->isHidden()) {
        mUi->actionEditJobQueue->setChecked(true);
    }

    this->show();
}

MainWindow::~MainWindow()
{
    SavePersistentSettings();
    delete mUi;
}

void MainWindow::LoadPersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("MainWindow");

    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    settings.endGroup();
}

void MainWindow::SavePersistentSettings()
{
    QSettings settings("siret", "molpher");
    settings.beginGroup("MainWindow");

    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());

    settings.endGroup();
}

void MainWindow::OnCloseBookmarks()
{
    mUi->actionBookmarks->setChecked(false);
}

void MainWindow::OnCloseJobsQueue()
{
    mUi->actionEditJobQueue->setChecked(false);
}

void MainWindow::OnCreateJob()
{
    CreateSnapshotDialog *snapshot = new CreateSnapshotDialog();
    snapshot->show();
}

void MainWindow::OnShowBookmarks(bool show)
{
    Bookmarks *b = gGlobalObjectsHolder.GetBookmarks();
    if (show) {
        b->show();
    } else {
        b->hide();
    }
}

void MainWindow::OnShowJobsQueue(bool show)
{
    if (show) {
        mJobsQueue->show();
    } else {
        mJobsQueue->hide();
    }
}

void MainWindow::OnCloseTab(int index)
{
    Tab *tab = qobject_cast<Tab *>(mUi->tabWidget->widget(index));
    mUi->tabWidget->removeTab(index);
    delete tab;
}

void MainWindow::OnOpenSnapshots()
{
    std::string storage("Temp/");
    storage += gFrontendId + "/";
    storage += boost::posix_time::to_iso_string(
        boost::posix_time::second_clock::local_time());

    try {
        boost::filesystem::create_directories(GenerateDirname(storage, -1));
    } catch (boost::filesystem::filesystem_error &exc) {
        SynchCout(exc.what());
    }

    std::vector<IterationSnapshotProxy> snapshots;
    QStringList filenames = GetFileNamesReadSNP(this);
    for (int i = 0; i < filenames.count(); ++i) {
        IterationSnapshot snp;
        if (ReadSnapshotFromFile(filenames.at(i).toStdString(), snp)) {
            WriteSnapshotToFile(GenerateFilename(storage, -1, i), snp);
            IterationSnapshotProxy proxy(storage, -1, i);
            snapshots.push_back(proxy);
        }
    }

    OnOpenAdhocTab(snapshots);
}

int MainWindow::GetLiveTabIndex(JobId jobId)
{
    Tab *tab = NULL;
    for (int i = 0; i < mUi->tabWidget->count(); ++i) {
        tab = qobject_cast<Tab *>(mUi->tabWidget->widget(i));
        if (tab->IsOnline() && (tab->GetJobId() == jobId)) {
            return i;
        }
    }
    return -1;
}

void MainWindow::OnOpenLiveTab(JobId jobId, IterationSnapshot &snapshot)
{
    int tabIndex = GetLiveTabIndex(jobId);

    if (-1 == tabIndex) {
        Tab *tab = new Tab(jobId, snapshot, false);
        tabIndex = mUi->tabWidget->addTab(tab, QString::number(jobId));

        connect(tab, SIGNAL(OpenAdhocTab(IterationSnapshotProxy &)),
            this, SLOT(OnOpenAdhocTab(IterationSnapshotProxy &)));
    }

    assert(-1 != tabIndex);

    mUi->tabWidget->setCurrentIndex(tabIndex);
}

void MainWindow::OnOpenDetachedTab(JobId jobId, IterationSnapshot &snapshot)
{
    QString label;
    label.append("\\");
    label.append(QString::number(snapshot.jobId));

    Tab *tab = new Tab(jobId, snapshot, true);
    int lastIndex = mUi->tabWidget->addTab(tab, label);
    mUi->tabWidget->setCurrentIndex(lastIndex);
}

void MainWindow::OnOpenAdhocTab(IterationSnapshotProxy &snapshot)
{
    QString label = QString::number(snapshot.jobId).append(" : ").append(
        QString::number(snapshot.iterIdx));

    std::vector<IterationSnapshotProxy> snapshots;
    snapshots.push_back(snapshot);
    Tab *tab = new Tab(snapshots);
    int lastIndex = mUi->tabWidget->addTab(tab, label);
    mUi->tabWidget->setCurrentIndex(lastIndex);
}

void MainWindow::OnOpenAdhocTab(std::vector<IterationSnapshotProxy> &snapshots)
{
    if (!snapshots.empty()) {
        QString label("~");
        Tab *tab = new Tab(snapshots);
        int lastIndex = mUi->tabWidget->addTab(tab, label);
        mUi->tabWidget->setCurrentIndex(lastIndex);
    }
}

void MainWindow::OnJobsUpdated(const JobGroup &jobs)
{
    Tab *tab = NULL;
    for (int i = 0; i < mUi->tabWidget->count(); ++i) {
        tab = qobject_cast<Tab *>(mUi->tabWidget->widget(i));
        bool finished = (jobs.mJobMap.find(tab->GetJobId()) == jobs.mJobMap.end());
        if (tab->IsOnline() && finished) {
            QString label = mUi->tabWidget->tabText(i);
            label.prepend("\\");
            mUi->tabWidget->setTabText(i, label);
            tab->OnDisplayOfflineState();
        }
    }
}

void MainWindow::OnUpdateJobHistory(JobId jobId, IterIdx maxIterIdx)
{
    Tab *tab = NULL;
    for (int i = 0; i < mUi->tabWidget->count(); ++i) {
        tab = qobject_cast<Tab *>(mUi->tabWidget->widget(i));
        if (tab->IsOnline() && (tab->GetJobId() == jobId)) {
            tab->FillSnapshotsList(0, maxIterIdx, true);
        }
    }
}

void MainWindow::OnConnect()
{
    ConnectDialog *connection = new ConnectDialog();
    connection->show();
}

void MainWindow::OnDisconnect()
{
    emit DisconnectFromBackend();
}

void MainWindow::OnAbout()
{
    QMessageBox::about(this, tr("About Molpher GUI"),
        tr("<h2>Molpher GUI</h2>"
        "<p>version <b>" MOLPH_VERSION "</b>"
        "<p>revision <b>" MOLPH_REVISION "</b> built on <b>" MOLPH_DATE ", " MOLPH_TIME "</b>"
        "<p>Copyright &copy; 2012 Vladimir Fiklik, Petr Koupy, Marek Mikes, Martin Straka, Peter Szepe"
        "<p>Copyright &copy; 2013 Petr Skoda"
        "<p>This program is free software: you can redistribute it and/or modify "
        "it under the terms of the GNU General Public License as published by "
        "the Free Software Foundation, either version 3 of the License, or "
        "(at your option) any later version."
        "<p>This program is distributed in the hope that it will be useful, "
        "but WITHOUT ANY WARRANTY; without even the implied warranty of "
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
        "GNU General Public License for more details."
        "<p>You should have received a copy of the GNU General Public License "
        "along with this program.  If not, see "
        "&lt;<a href=\"http://www.gnu.org/licenses/\">http://www.gnu.org/licenses/</a>&gt;."
        "<h4>Acknowledgements</h4>"
        "<p>This program was developed under supervision of:<br>"
        "RNDr. David Hoksza, Ph.D. | SIRET research group, Faculty of Mathematics and Physics, Charles University in Prague<br>"
        "Doc. Daniel Svozil, Ph.D. | Laboratory of Informatics and Chemistry, Institute of Chemical Technology, Prague"
        "<p>This program uses following GNU GPL compatible third-party software:<br>"
        "Qt by Nokia Corporation<br>"
        "Boost by boost.org<br>"
        "RCF by Delta V Software<br>"
        "RDKit by Rational Discovery LLC<br>"
        "indigo-depict by GGA Software Services LLC<br>"
        "zlib by Jean-loup Gailly and Mark Adler"
        ));
}

void MainWindow::OnLegend()
{
    QDialog *legendDialog = new QDialog(this);
    legendDialog->setWindowTitle(tr("Legend"));
    legendDialog->setWindowFlags(
        (legendDialog->windowFlags() & ~Qt::WindowContextHelpButtonHint) |
        Qt::WindowStaysOnTopHint);
    legendDialog->setAttribute(Qt::WA_DeleteOnClose, true);
    QLabel *label = new QLabel(legendDialog);
    QPixmap pixmap;
    if (pixmap.load(":/images/legend.png")) {
        label->setPixmap(pixmap);
    } else {
        label->setText(tr("\n\tLegend is unavailable.\t\n"));
    }

    QGridLayout *layout = new QGridLayout(legendDialog);
    layout->addWidget(label);
    layout->setMargin(0);
    legendDialog->setLayout(layout);
    legendDialog->layout()->setSizeConstraint(
        QLayout::SetFixedSize);

    legendDialog->show();
}

void MainWindow::OnDisplayOnlineState(const QString &ip)
{
    mConnectionState.setText(tr("Online (").append(ip).append(")"));
    mUi->actionConnect->setEnabled(false);
    mUi->actionDisconnect->setEnabled(true);
    mUi->actionNewJob->setEnabled(true);
}

void MainWindow::OnDisplayOfflineState()
{
    for (int i = 0; i < mUi->tabWidget->count(); ++i) {
        QString label = mUi->tabWidget->tabText(i);
        bool wasOnline = !label.contains("\\") &&
            !label.contains(":") && !label.contains("~");
        if (wasOnline) {
            label.prepend("\\");
            mUi->tabWidget->setTabText(i, label);
        }
    }

    mConnectionState.setText(tr("Offline"));
    mUi->actionConnect->setEnabled(true);
    mUi->actionDisconnect->setEnabled(false);
    mUi->actionNewJob->setEnabled(false);
    JobGroup empty;
    mJobsQueue->OnUpdateJobs(empty);
}

void MainWindow::OnClearSettings() {
    // delete all stored settings
    QSettings settings("siret", "molpher");
    settings.clear();
    // write the result on disk
    settings.sync();
}
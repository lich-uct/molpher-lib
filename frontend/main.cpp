/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Vladimir Fiklik

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

#include <iostream>
#include <string>

#include <RCF/RCF.hpp>

#include <boost/program_options.hpp>

#include <QtCore/QObject>
#include <QtGui/QApplication>
#include <QtGui/QStatusBar>

#include "inout.h"
#include "FrontendCommunicator.h"
#include "form/MainWindow.h"
#include "tests/ChangeSelectorsTester.h"
#include "tests/DatabaseTester.h"
#include "tests/DecoysTester.h"
#include "tests/NeighborhoodTester.h"
#include "tests/SleepWakeTester.h"
#include "tests/PruningTester.h"
#include "tests/BasicTester.h"
#include "tests/TestManager.h"

void RegisterFastTesters(TestManager &manager)
{
    ChangeSelectorsTester *changeSelTester =  new ChangeSelectorsTester(false, "pass");
    manager.RegisterTester(changeSelTester);

    SleepWakeTester *sleepWakeTester = new SleepWakeTester("pass2");
    manager.RegisterTester(sleepWakeTester);
}

void RegisterSlowTesters(TestManager &manager)
{
    //DatabaseTester *databaseTester = new DatabaseTester();
    //manager.RegisterTester(databaseTester);

    PruningTester *pruningTester = new PruningTester("pass");
    manager.RegisterTester(pruningTester);

    ChangeSelectorsTester *changeSelTester = new ChangeSelectorsTester(true, "pass3");
    manager.RegisterTester(changeSelTester);

    //DecoysTester *decoysTester = new DecoysTester("pass4");
    //manager.RegisterTester(decoysTester);

    NeighborhoodTester *neighbTester = new NeighborhoodTester();
    manager.RegisterTester(neighbTester);
}

void RegisterBasicTester(TestManager &manager)
{
    BasicTester *basicTester = new BasicTester();
    manager.RegisterTester(basicTester);
}

void Run(int argc, char *argv[])
{
    std::cout << "Initializing..." << std::endl;

    QApplication app(argc, argv);

    if (gCommunicator.Init()) {
        MainWindow mainWindow;
        QObject::connect(&gCommunicator, SIGNAL(ReportStatus(const QString &, int)),
            mainWindow.statusBar(), SLOT(showMessage(const QString &, int)));
        gCommunicator.ConnectToBackend("localhost"); // Try localhost.
        SynchCout(std::string("Frontend initialized."));
        app.exec();

        gCommunicator.DisconnectFromBackend();
        gCommunicator.DeInit();
    } else {
        std::cout << "Cannot initialize communicator." << std::endl;
    }

    std::cout << "Frontend terminated." << std::endl;
}

void Test(int argc, char *argv[], char testOption)
{
    std::cout << "Initializing test mode..." << std::endl;

    QApplication app(argc, argv);

    if (gCommunicator.Init()) {

        TestManager testManager;

        switch (testOption) {
        case 's':
            RegisterSlowTesters(testManager);
            break;
        case 'f':
            RegisterFastTesters(testManager);
            break;
        case 'b':
            RegisterBasicTester(testManager);
            break;
        case 'a':
            RegisterFastTesters(testManager);
            RegisterSlowTesters(testManager);
            RegisterBasicTester(testManager);
            break;
        default:
            break;
        }

        gCommunicator.ConnectToBackend("localhost");
        SynchCout(std::string("Test mode initialized."));
        SynchCout(std::string("Note: Backend must be running."));
        app.exec();

        gCommunicator.DisconnectFromBackend();
        gCommunicator.DeInit();
    } else {
        std::cout << "Cannot initialize communicator." << std::endl;
    }

    std::cout << "Frontend terminated." << std::endl;
}

int main(int argc, char *argv[])
{
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Show help")
        ("test,T",boost::program_options::value<char>(), "Storage path for results")
            ;

    boost::program_options::variables_map varMap;
    try {
        boost::program_options::store(
            boost::program_options::parse_command_line(argc, argv, desc), varMap);
    } catch (boost::program_options::error &exc) {
        std::cout << desc << std::endl;
        return -1;
    }
    boost::program_options::notify(varMap);

    if (varMap.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    bool test = false;

    char testOption;
    if (varMap.count("test")) {
        test = true;
        testOption = varMap["test"].as<char>();
    }

    RCF::init();
    if (test) {
        Test(argc, argv, testOption);
    } else {
        Run(argc, argv);
    }
    RCF::deinit();

    return 0;
}

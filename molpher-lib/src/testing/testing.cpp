
#include <string>
#include <iostream>
#include <thread>

#include <tr1/functional>
#include <tbb/task.h>
#include <tbb/compat/thread>

#include "IterationSnapshot.h"
#include "core/JobManager.h"
#include "path_finders/BasicPathFinder.hpp"
#include "inout.h"
#include "extensions/SAScore.h"

#include "testing/testing.hpp"

void run_path_finder(
    std::string &storagePath
    , std::string &jobFile
    , int threadCnt
) {
    
    SAScore::loadData(); // load data for prediction of synthetic feasibility
    
    tbb::task_group_context pathFinderTbbCtx;
    std::string dummy;
    JobManager jobManager(&pathFinderTbbCtx, storagePath, dummy, false);
    jobManager.AddJobFromFile(jobFile);
    BasicPathFinder pathFinder(&pathFinderTbbCtx, &jobManager, threadCnt);
    std::thread pathFinderThread(std::tr1::ref(pathFinder));
    SynchCout(std::string("Backend initialized.\nWorking..."));
    pathFinderThread.join();
    SynchCout(std::string("Halting..."));
    jobManager.Halt();
    std::cout << "Backend terminated." << std::endl;
    
    SAScore::destroyInstance(); // should free data, maybe not necessary
}
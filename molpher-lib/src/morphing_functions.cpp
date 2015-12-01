
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

#include "morphing_functions.hpp"


void run_path_finder(
    const std::string &storagePath
    , const std::string &jobFile
    , int threadCnt
) {
    
    SAScore::loadData();
    
    tbb::task_group_context pathFinderTbbCtx;
    std::string dummy;
    std::string storagePathnonconst(storagePath);
    JobManager jobManager(&pathFinderTbbCtx, storagePathnonconst, dummy, false);
    std::string jobfilenonconst(jobFile);
    jobManager.AddJobFromFile(jobfilenonconst);
    BasicPathFinder pathFinder(&pathFinderTbbCtx, &jobManager, threadCnt);
    std::thread pathFinderThread(std::tr1::ref(pathFinder));
    SynchCout(std::string("Backend initialized.\nWorking..."));
    pathFinderThread.join();
    SynchCout(std::string("Halting..."));
    jobManager.Halt();
    std::cout << "Backend terminated." << std::endl;
}



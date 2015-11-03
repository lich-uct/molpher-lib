
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

/** 
 * Runs the original morphing algorithm 
 * using the BasicPathFinder implementation.
 * 
 * Note: Don't forget to call SAScore::loadData() before the call!
 * 
 * @param storagePath name of the folder where the results are saved
 * @param jobFile path to a job file or XML template
 * @param thredCnt max. number of threads to create
 */
void run_path_finder(
    const std::string &storagePath
    , const std::string &jobFile
    , int threadCnt
) {
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



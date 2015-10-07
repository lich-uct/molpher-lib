
#include <stdexcept>

#include "MorphingManager.hpp"
#include "inout.h"

//MoprhingManager::MoprhingManager() :
//    mSnapshot() 
//{
//    initContext();
//}

MoprhingManager::MoprhingManager(const std::string &jobFile) :
MoprhingManager(jobFile, -1) {
    // no action
}

MoprhingManager::MoprhingManager(const std::string &jobFile, int threadCnt) :
mThreadCount(threadCnt) {
    if (ReadSnapshotFromFile(jobFile, mSnapshot) && mSnapshot.IsValid()) {
        SynchCout(std::string("Snapshot successfully created from: ").append(jobFile));
    } else {
        SynchCout(std::string("Cannot open job: ").append(jobFile));
        throw std::runtime_error("Job file cannot be opened.");
    }
    initContext();
}

void MoprhingManager::initContext() {
    PathFinderContext::SnapshotToContext(mSnapshot, mContext);
}

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "FindLeavesOper.hpp"

FindLeavesOper::FindLeaves::FindLeaves(ExplorationTree::MoleculeVector &leaves) :
mLeaves(leaves) {
}

void FindLeavesOper::FindLeaves::operator()(
        const PathFinderContext::CandidateMap::range_type &candidates) const {
    PathFinderContext::CandidateMap::iterator it;
    for (it = candidates.begin(); it != candidates.end(); it++) {
        if (!it->second.parentSmile.empty()) {
            it->second.itersWithoutDistImprovement++;
        }
        bool isLeaf = it->second.descendants.empty();
        if (isLeaf) {
            mLeaves.push_back(it->second);
        }
    }
}

FindLeavesOper::FindLeavesOper(ExplorationTree& expTree) : TreeOperation(expTree) {
    // no action
}

void FindLeavesOper::operator()() {
    PathFinderContext& ctx = TreeOperation::fetchTreeContext();
    tbb::task_group_context tbbCtx;
    
    tbb::task_scheduler_init scheduler;
    scheduler.terminate();
    scheduler.initialize(2); // TODO expose as parameter
    
    FindLeaves findLeaves(leaves);
    tbb::parallel_for(
                PathFinderContext::CandidateMap::range_type(ctx.candidates),
                findLeaves, tbb::auto_partitioner(), tbbCtx);
//        stageStopwatch.ReportElapsedMiliseconds("FindLeaves", true);
}

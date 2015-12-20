
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "molpher_API/operations/FindLeavesOper.hpp"

FindLeavesOper::FindLeaves::FindLeaves(ExplorationTree::MoleculePointerVector &leaves, bool increment_iters_without_dist_improve) :
mLeaves(leaves) 
, mIncrementDistImproveCounter(increment_iters_without_dist_improve)
{
}

void FindLeavesOper::FindLeaves::operator()(
        const PathFinderContext::CandidateMap::range_type &candidates) const {
    PathFinderContext::CandidateMap::iterator it;
    for (it = candidates.begin(); it != candidates.end(); it++) {
        if (!it->second.parentSmile.empty() && mIncrementDistImproveCounter) {
            it->second.itersWithoutDistImprovement++;
        }
        bool isLeaf = it->second.descendants.empty();
        if (isLeaf) {
            mLeaves.push_back(&(it->second));
        }
    }
}

FindLeavesOper::FindLeavesOper(ExplorationTree& expTree, bool increment_iters_without_dist_improve) : 
TreeOperation(expTree)
, mIncrementDistImproveCounter(increment_iters_without_dist_improve) 
{
    // no action
}

FindLeavesOper::FindLeavesOper() : 
TreeOperation()
, mIncrementDistImproveCounter(false) 
{
    // no action
}

void FindLeavesOper::operator()() {
    PathFinderContext& ctx = TreeOperation::fetchTreeContext();
    tbb::task_group_context tbbCtx;
    
    tbb::task_scheduler_init scheduler;
    if (threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(threadCnt);
    }
    
    FindLeaves findLeaves(leaves, mIncrementDistImproveCounter);
    tbb::parallel_for(
                PathFinderContext::CandidateMap::range_type(ctx.candidates),
                findLeaves, tbb::auto_partitioner(), tbbCtx);
//        stageStopwatch.ReportElapsedMiliseconds("FindLeaves", true);
}

const std::vector<MolpherMol>& FindLeavesOper::fetchLeaves() {
    if (this->tree) {
        std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
        for (auto leaf : leaves) {
            ret->push_back(MolpherMol(*leaf));
        }
        return *ret;
    } else {
        throw std::runtime_error("Cannot find leaves. No tree associated with this instance.");
    }
}

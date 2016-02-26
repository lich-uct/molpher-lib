
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "operations/FindLeavesOper.hpp"
#include "FindLeavesOperImpl.hpp"
#include "tbb/concurrent_hash_map.h"
#include "core/API/MolpherMolImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "core/misc/global_types.h"

FindLeavesOper::FindLeavesOper(bool increment_iters_without_dist_improve) : 
pimpl(new FindLeavesOper::FindLeavesOperImpl(increment_iters_without_dist_improve))
{
    setTreeOperPimpl(pimpl);
}

FindLeavesOper::FindLeavesOper(std::shared_ptr<ExplorationTree> expTree, bool increment_iters_without_dist_improve) :
pimpl(new FindLeavesOper::FindLeavesOperImpl(expTree, increment_iters_without_dist_improve))
{
    setTreeOperPimpl(pimpl);
}

void FindLeavesOper::operator()() {
    (*pimpl)();
}

std::vector<std::shared_ptr<MolpherMol> > FindLeavesOper::fetchLeaves() {
    return *(pimpl->fetchLeaves());
}

std::unique_ptr<MolVector> FindLeavesOper::FindLeavesOperImpl::fetchLeaves() {
    auto ret = std::unique_ptr<MolVector>(new MolVector());
    if (!getTree()) {
        std::cerr << "WARNING: No tree specified. No leaves acquired" << std::endl;
        return ret;
    } else {
        concurrent_vector_to_vector<std::shared_ptr<MolpherMol> >(leaves, *ret);
        return ret;
    }
}

FindLeavesOper::FindLeavesOperImpl::FindLeavesOperImpl(
    std::shared_ptr<ExplorationTree> expTree
    , bool increment_iters_without_dist_improve) 
    :
    TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
    , incrementDistImproveCounter(increment_iters_without_dist_improve)
{
    // no action
}

FindLeavesOper::FindLeavesOperImpl::FindLeavesOperImpl(bool increment_iters_without_dist_improve) : 
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, incrementDistImproveCounter(increment_iters_without_dist_improve)
{
    // no action
}


FindLeavesOper::FindLeavesOperImpl::FindLeaves::FindLeaves(ConcurrentMolVector &leaves, bool increment_iters_without_dist_improve) :
mLeaves(leaves) 
, mIncrementDistImproveCounter(increment_iters_without_dist_improve)
{
}

void FindLeavesOper::FindLeavesOperImpl::FindLeaves::operator()(
        const TreeMap::range_type &tree) const {
    TreeMap::iterator it;
    
    for (it = tree.begin(); it != tree.end(); it++) {
        if (!it->second->getParentSMILES().empty() && mIncrementDistImproveCounter) {
            it->second->increaseItersWithoutDistImprovement();
        }
        bool isLeaf = it->second->getDescendants().empty();
        if (isLeaf) {
            mLeaves.push_back(it->second);
        }
    }
}

void FindLeavesOper::FindLeavesOperImpl::operator()() {
    tbb::task_group_context tbbCtx;
    
    tbb::task_scheduler_init scheduler;
    auto tree_pimpl = getTree()->pimpl;
    if (tree_pimpl->threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(tree_pimpl->threadCnt);
    }
    
    FindLeavesOper::FindLeavesOperImpl::FindLeaves findLeaves(leaves, incrementDistImproveCounter);
    tbb::parallel_for(
                TreeMap::range_type(tree_pimpl->treeMap),
                findLeaves, tbb::auto_partitioner(), tbbCtx);
}

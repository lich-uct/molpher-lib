
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "operations/FindLeavesOper.hpp"
#include "FindLeavesOperImpl.hpp"
#include "tbb/concurrent_hash_map.h"
#include "core/API/MolpherMolImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"

FindLeavesOper::FindLeavesOper() : 
TreeOperation() 
, pimpl(new FindLeavesOper::FindLeavesOperImpl())
{
    // no action
}

FindLeavesOper::FindLeavesOper(std::shared_ptr<ExplorationTree> expTree, bool increment_iters_without_dist_improve = false) :
TreeOperation(expTree)
, pimpl(new FindLeavesOper::FindLeavesOperImpl(expTree->pimpl, increment_iters_without_dist_improve))
{
    // no action
}

void FindLeavesOper::operator()() {
    (*pimpl)();
}

std::vector<std::shared_ptr<MolpherMol> > FindLeavesOper::fetchLeaves() {
    return pimpl->fetchLeaves();
}

std::shared_ptr<MolVectorAPI> FindLeavesOper::FindLeavesOperImpl::fetchLeaves() {
    std::shared_ptr<MolVectorAPI> ret;
    if (!getTree()) {
        std::cerr << "WARNING: No tree specified. No leaves aquired" << std::endl;
        return ret;
    } else {
        MolVector ret_pimpl;
        concurrent_vector_to_vector<std::shared_ptr<MolpherMol::MolpherMolImpl>>(leaves, ret_pimpl);
        for (auto& mol_pimpl : ret_pimpl) {
            ret->push_back(std::make_shared<MolpherMol>(mol_pimpl));
        }
        return ret;
    }
}

FindLeavesOper::FindLeavesOperImpl::FindLeavesOperImpl(
    std::shared_ptr<ExplorationTree::ExplorationTreeImpl> expTree
    , bool increment_iters_without_dist_improve) 
    :
    TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
    , mIncrementDistImproveCounter(increment_iters_without_dist_improve)
{
    // no action
}

FindLeavesOper::FindLeavesOperImpl::FindLeavesOperImpl() : 
mIncrementDistImproveCounter(false)
{
    // no action
}


FindLeavesOper::FindLeavesOperImpl::FindLeaves::FindLeaves(ConcurrentMolVector &leaves, bool increment_iters_without_dist_improve) :
mLeaves(leaves) 
, mIncrementDistImproveCounter(increment_iters_without_dist_improve)
{
}

void FindLeavesOper::FindLeavesOperImpl::FindLeaves::operator()(
        const TreeMap &tree) const {
    TreeMap::iterator it;
    
    for (it = tree.begin(); it != tree.end(); it++) {
        if (!it->second->getParentSMILES().empty() && mIncrementDistImproveCounter) {
            it->second->getItersWithoutDistImprovement()++;
        }
        bool isLeaf = it->second->getDescendants()->empty();
        if (isLeaf) {
            mLeaves.push_back(&(it->second));
        }
    }
}

void FindLeavesOper::FindLeavesOperImpl::operator()() {
    auto& tree = getTree();
    tbb::task_group_context tbbCtx;
    
    tbb::task_scheduler_init scheduler;
    if (tree->threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(tree->threadCnt);
    }
    
    FindLeaves findLeaves(leaves, mIncrementDistImproveCounter);
    tbb::parallel_for(
                TreeMap::range_type(tree->treeMap),
                findLeaves, tbb::auto_partitioner(), tbbCtx);
}

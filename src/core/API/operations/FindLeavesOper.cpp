
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "operations/FindLeavesOper.hpp"
#include "FindLeavesOperImpl.hpp"

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
    MolVector ret_pimpl;
    concurrent_vector_to_vector<std::shared_ptr<MolpherMol::MolpherMolImpl>>(leaves, ret_pimpl);
    MolVectorAPI ret;
    for (auto& mol_pimpl : ret_pimpl) {
        ret->push_back(std::make_shared<MolpherMol>(mol_pimpl));
    }
    return ret;
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

FindLeavesOper::FindLeavesOperImpl::FindLeavesOperImpl() {
    // TODO: next
}

void FindLeavesOper::FindLeavesOperImpl::operator()() {
    // TODO: next
}










//FindLeavesOper::FindLeaves::FindLeaves(ExplorationTree::MoleculePointerVector &leaves, bool increment_iters_without_dist_improve) :
//mLeaves(leaves) 
//, mIncrementDistImproveCounter(increment_iters_without_dist_improve)
//{
//}
//
//void FindLeavesOper::FindLeaves::operator()(
//        const PathFinderContext::CandidateMap::range_type &candidates) const {
//    PathFinderContext::CandidateMap::iterator it;
//    for (it = candidates.begin(); it != candidates.end(); it++) {
//        if (!it->second.parentSmile.empty() && mIncrementDistImproveCounter) {
//            it->second.itersWithoutDistImprovement++;
//        }
//        bool isLeaf = it->second.descendants.empty();
//        if (isLeaf) {
//            mLeaves.push_back(&(it->second));
//        }
//    }
//}
//
//FindLeavesOper::FindLeavesOper(ExplorationTree& expTree, bool increment_iters_without_dist_improve) : 
//TreeOperation(expTree)
//, mIncrementDistImproveCounter(increment_iters_without_dist_improve) 
//{
//    // no action
//}
//
//FindLeavesOper::FindLeavesOper() : 
//TreeOperation()
//, mIncrementDistImproveCounter(false) 
//{
//    // no action
//}
//
//void FindLeavesOper::operator()() {
//    PathFinderContext& ctx = TreeOperation::fetchTreeContext();
//    tbb::task_group_context tbbCtx;
//    
//    tbb::task_scheduler_init scheduler;
//    if (threadCnt > 0) {
//        scheduler.terminate();
//        scheduler.initialize(threadCnt);
//    }
//    
//    FindLeaves findLeaves(leaves, mIncrementDistImproveCounter);
//    tbb::parallel_for(
//                PathFinderContext::CandidateMap::range_type(ctx.candidates),
//                findLeaves, tbb::auto_partitioner(), tbbCtx);
////        stageStopwatch.ReportElapsedMiliseconds("FindLeaves", true);
//}
//
//const std::vector<MolpherMol>& FindLeavesOper::fetchLeaves() {
//    if (this->tree) {
//        std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
//        for (auto leaf : leaves) {
//            ret->push_back(MolpherMol(*leaf));
//        }
//        return *ret;
//    } else {
//        throw std::runtime_error("Cannot find leaves. No tree associated with this instance.");
//    }
//}

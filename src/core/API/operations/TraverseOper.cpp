
#include "operations/TraverseOper.hpp"
#include "TraverseOperImpl.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "TreeOperationImpl.hpp"

TraverseOper::TraverseOper(std::shared_ptr<ExplorationTree> expTree, std::shared_ptr<TraverseCallback> callback) : 
pimpl(new TraverseOper::TraverseOperImpl(expTree, callback)) 
{
    setTreeOperPimpl(pimpl);
}

TraverseOper::TraverseOper(std::shared_ptr<TraverseCallback> callback) : 
pimpl(new TraverseOper::TraverseOperImpl(callback))
{
    setTreeOperPimpl(pimpl);
}

TraverseOper::TraverseOper(
    std::shared_ptr<ExplorationTree> expTree
    , std::shared_ptr<TraverseCallback> callback
    , const std::string& rootSMILES
    ) 
: 
pimpl(new TraverseOper::TraverseOperImpl(expTree, callback, rootSMILES))  
{
    setTreeOperPimpl(pimpl);
}

void TraverseOper::operator()() {
    (*pimpl)();
}


// pimpl

TraverseOper::TraverseOperImpl::TraverseOperImpl(
    std::shared_ptr<ExplorationTree> expTree
    , std::shared_ptr<TraverseCallback> callback
    , const std::string& rootSMILES
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, callback(callback)
, root(expTree->fetchMol(rootSMILES))
{
    // no action
}

TraverseOper::TraverseOperImpl::TraverseOperImpl(
    std::shared_ptr<ExplorationTree> expTree
    , std::shared_ptr<TraverseCallback> callback
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, callback(callback)
{
    // no action
}

TraverseOper::TraverseOperImpl::TraverseOperImpl(std::shared_ptr<TraverseCallback> callback) :
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, callback(callback)
{
    // no action
}

TraverseOper::TraverseOperImpl::TraversalFunctor::TraversalFunctor(
    std::shared_ptr<ExplorationTree> tree
    , std::shared_ptr<TraverseCallback> callback
    ) 
: 
mTree(tree)
, mCallback(callback)
{
    // no action
}

void TraverseOper::TraverseOperImpl::TraversalFunctor::operator()(const std::string& smile, tbb::parallel_do_feeder<std::string>& feeder) const {
    auto tree_pimpl = mTree->pimpl;

    TreeMap::accessor ac;
    tree_pimpl->treeMap.find(ac, smile);
    assert(!ac.empty());

    makeCallback(ac->second);

    std::set<std::string>::const_iterator it;
    auto descendants = ac->second->getDescendants();
    for (it = descendants.begin();
            it != descendants.end(); it++) {
        feeder.add(*it);
    }
}

void TraverseOper::TraverseOperImpl::TraversalFunctor::makeCallback(std::shared_ptr<MolpherMol> morph) const {
    (*mCallback)(morph);
}

void TraverseOper::TraverseOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_pimpl = tree->pimpl;
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (tree_pimpl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree_pimpl->threadCnt);
        }
        ConcurrentSmileVector queue;
        if (!root) {
            root = tree->fetchMol(tree_pimpl->source.getSMILES());
        }
        queue.push_back(root->getSMILES());

        TraversalFunctor functor(tree, callback);
        tbb::parallel_do(queue.begin(), queue.end(), functor, tbbCtx);
    } else {
        throw std::runtime_error("Cannot traverse the tree. None associated with this instance.");
    }
}



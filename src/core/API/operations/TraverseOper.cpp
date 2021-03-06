/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/misc/inout.h"

#include "operations/TraverseOper.hpp"
#include "TraverseOperImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"

TraverseOper::TraverseOper(std::shared_ptr<ExplorationTree> expTree, TraverseCallback& callback) : 
pimpl(new TraverseOper::TraverseOperImpl(expTree, callback)) 
{
    setTreeOperPimpl(pimpl);
}

TraverseOper::TraverseOper(TraverseCallback& callback) : 
pimpl(new TraverseOper::TraverseOperImpl(callback))
{
    setTreeOperPimpl(pimpl);
}

TraverseOper::TraverseOper(
    std::shared_ptr<ExplorationTree> expTree
    , TraverseCallback& callback
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
    , TraverseCallback& callback
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
    , TraverseCallback& callback
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, callback(callback)
{
    // no action
}

TraverseOper::TraverseOperImpl::TraverseOperImpl(TraverseCallback& callback) :
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, callback(callback)
{
    // no action
}

TraverseOper::TraverseOperImpl::TraversalFunctor::TraversalFunctor(
    std::shared_ptr<ExplorationTree> tree
    , TraverseCallback& callback
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

    std::shared_ptr<MolpherMol> mol = ac->second;
//    int x = mol.use_count();
    ac.release();
    makeCallback(mol);
//    SynchCout("use count (before/after): " + parseNumber(x) + "/" + parseNumber(mol.use_count()));

    std::set<std::string>::const_iterator it;
    auto descendants = mol->getDescendants();
    for (it = descendants.begin();
            it != descendants.end(); it++) {
        feeder.add(*it);
    }
}

void TraverseOper::TraverseOperImpl::TraversalFunctor::makeCallback(std::shared_ptr<MolpherMol> morph) const {
    mCallback(morph);
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
            root = tree->fetchMol(tree_pimpl->source->getSMILES());
        }
        queue.push_back(root->getSMILES());

        TraversalFunctor functor(tree, callback);
        tbb::parallel_do(queue.begin(), queue.end(), functor, tbbCtx);
    } else {
        throw std::runtime_error("Cannot traverse the tree. None associated with this instance.");
    }
}



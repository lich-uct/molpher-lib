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

#include "operations/PruneTreeOper.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "PruneTreeOperImpl.hpp"
#include "TreeOperationImpl.hpp"

PruneTreeOper::PruneTreeOper(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new PruneTreeOper::PruneTreeOperImpl(expTree))
{
    setTreeOperPimpl(pimpl);
}

PruneTreeOper::PruneTreeOper() :
pimpl(new PruneTreeOper::PruneTreeOperImpl())
{
    setTreeOperPimpl(pimpl);
}

void PruneTreeOper::operator()() {
    (*pimpl)();
}

// pimpl

PruneTreeOper::PruneTreeOperImpl::PruneTreeOperImpl(std::shared_ptr<ExplorationTree> expTree) :
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
{
    // no action
}

PruneTreeOper::PruneTreeOperImpl::PruneTreeOperImpl() :
TreeOperation::TreeOperationImpl::TreeOperationImpl()
{
    // no action
}



PruneTreeOper::PruneTreeOperImpl::PruneTree::PruneTree(std::shared_ptr<ExplorationTree> expTree) :
mTree(expTree)
{
    // no action
}

void PruneTreeOper::PruneTreeOperImpl::PruneTree::operator()(const std::string& smile, tbb::parallel_do_feeder<std::string>& feeder) const {
    TreeMap::accessor ac;
    auto tree_pimpl = mTree->pimpl;
    tree_pimpl->treeMap.find(ac, smile);
    assert(!ac.empty());

    bool prune = (ac->second->getItersWithoutDistImprovement() > tree_pimpl->params.itThreshold);
    if (prune) {

        bool tooManyDerivations = false;
        MorphDerivationMap::const_accessor acDerivations;
        if (tree_pimpl->morphDerivations.find(acDerivations, smile)) {
            tooManyDerivations = (acDerivations->second > tree_pimpl->params.cntMaxMorphs);
        }
        
        acDerivations.release();
        ac.release();
        
        if (tooManyDerivations) {
            eraseSubtree(smile, false);
            SynchCout("Pruned: " + smile);
        } else {
            eraseSubtree(smile, true);
            SynchCout("Pruned (descendents only): " + smile);
        }

    } else {
        std::set<std::string>::const_iterator it;
        auto descendants = ac->second->getDescendants();
        for (it = descendants.begin();
                it != descendants.end(); it++) {
            feeder.add(*it);
        }
    }
}

void PruneTreeOper::PruneTreeOperImpl::PruneTree::eraseSubtree(const std::string& smile, bool descendents_only) const {
    mTree->deleteSubtree(smile, descendents_only);
}


void PruneTreeOper::PruneTreeOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (tree->pimpl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree->pimpl->threadCnt);
        }
        ConcurrentSmileVector queue;
        queue.push_back(tree->pimpl->source.getSMILES());

        PruneTree functor(tree);
        tbb::parallel_do(queue.begin(), queue.end(), functor, tbbCtx);
    } else {
        throw std::runtime_error("Cannot prune. No tree attached to this instance.");
    }
}

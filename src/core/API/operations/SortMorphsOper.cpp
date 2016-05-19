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

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_sort.h>

#include "operations/SortMorphsOper.hpp"
#include "SortMorphsOperImpl.hpp"
#include "TreeOperationImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"

SortMorphsOper::SortMorphsOper(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new SortMorphsOper::SortMorphsOperImpl(expTree))
{
     setTreeOperPimpl(pimpl);
}

SortMorphsOper::SortMorphsOper() :
pimpl(new SortMorphsOper::SortMorphsOperImpl())
{
    setTreeOperPimpl(pimpl);
}

void SortMorphsOper::operator()() {
    (*pimpl)();
}

SortMorphsOper::SortMorphsOperImpl::SortMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree) :
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
{
    // no action
}

SortMorphsOper::SortMorphsOperImpl::SortMorphsOperImpl() :
TreeOperation::TreeOperationImpl::TreeOperationImpl()
{
    // no action
}

bool SortMorphsOper::SortMorphsOperImpl::CompareMorphs::operator()(
        const std::shared_ptr<MolpherMol> &a, const std::shared_ptr<MolpherMol> &b) const {
    /* Morphs are rated according to their proximity to the connecting line
     between their closest decoy and the target (i.e. sum of both distances is
     minimal on the connecting line between decoy and target). When sums for
     both morphs are equal, it is possible (but not necessary) that both
     morphs lie on the same connecting line. In that case, morphs are
     rated only according to their proximity to the target. Such comparison
     should allow convergence to the target even in the late stages of the
     algorithm when majority of morphs lie on the connecting line between
     decoy closest to the target and the target itself. */

    double aSum = a->getDistToTarget();
    double bSum = b->getDistToTarget();

    bool approximatelyEqual = (
            fabs(aSum - bSum) <= (32 * DBL_EPSILON * fmax(fabs(aSum), fabs(bSum))));

    if (approximatelyEqual) {
        return a->getDistToTarget() < b->getDistToTarget();
    } else {
        return aSum < bSum;
    }

    /**
     * Just sort based on distance to decoy .. or target. We prefer those
     * with target.
     */

    /* Experimental evaluation
    if (a.nextDecoy == -1 && b.nextDecoy == -1) {
        // both go for target, so take just their distance
        return a.distToTarget < b.distToTarget;
    } else if (a.nextDecoy == -1) {
        // a goes for target while b not -> a is greater
        return true;
    } else if (b.nextDecoy == -1) {
        // b goes for target while a not -> b is greater
        return false;
    } else {
        // decide based on nextDecoy
        if (a.nextDecoy == b.nextDecoy) {
            // go for same decoy use distance to closes decoy
            return a.distToClosestDecoy < b.distToClosestDecoy;
        } else {
            // greater is better
            return a.nextDecoy > b.nextDecoy;
        }
    }*/
}

void SortMorphsOper::SortMorphsOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_pimpl = tree->pimpl;
        tbb::task_scheduler_init scheduler;
        if (tree_pimpl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree_pimpl->threadCnt);
        }

        CompareMorphs compareMorphs;
        tbb::parallel_sort(tree_pimpl->candidates.begin(), tree_pimpl->candidates.end(), compareMorphs);
    } else {
        throw std::runtime_error("Cannot sort morphs. No tree attached to this instance.");
    }
}


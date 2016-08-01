/*
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
#include <tbb/parallel_for.h>

#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"

#include "operations/CleanMorphsOper.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "CleanMorphsOperImpl.hpp"

CleanMorphsOper::CleanMorphsOper(std::shared_ptr<ExplorationTree> expTree) :
		pimpl(new CleanMorphsOper::CleanMorphsOperImpl(expTree))
{
	setTreeOperPimpl(pimpl);
}

CleanMorphsOper::CleanMorphsOper() :
		pimpl(new CleanMorphsOper::CleanMorphsOperImpl())
{
	setTreeOperPimpl(pimpl);
}

void CleanMorphsOper::operator()() {
	(*pimpl)();
}

// pimpl

CleanMorphsOper::CleanMorphsOperImpl::CleanMorphsOperImpl(
		std::shared_ptr<ExplorationTree> expTree
)
		:
		TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
{
	// no action
}

CleanMorphsOper::CleanMorphsOperImpl::CleanMorphsOperImpl()
		:
		TreeOperation::TreeOperationImpl::TreeOperationImpl()
{
	// no action
}

CleanMorphsOper::CleanMorphsOperImpl::CleanMorphs::CleanMorphs(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_impl,
																   size_t globalMorphCount, ConcurrentMolVector &morphs, std::vector<bool> &survivors, ConcurrentMolVector &new_candidates
) :
		mTreePimpl(tree_impl),
		mGlobalMorphCount(globalMorphCount),
		mMorphs(morphs),
		mSurvivors(survivors),
		mNewCandidates(new_candidates)
{
	assert(mMorphs.size() == mSurvivors.size());
}

void CleanMorphsOper::CleanMorphsOperImpl::CleanMorphs::operator()(const tbb::blocked_range<size_t> &r) const {

	for (size_t idx = r.begin(); idx != r.end(); ++idx) {

		if (mSurvivors[idx]) {
			mNewCandidates.push_back(mMorphs[idx]);
		}
	}
}

void CleanMorphsOper::CleanMorphsOperImpl::operator()() {
	auto tree = getTree();
	if (tree) {
		auto tree_impl = tree->pimpl;
		tbb::task_group_context tbbCtx;
		tbb::task_scheduler_init scheduler;
		if (tree_impl->threadCnt > 0) {
			scheduler.terminate();
			scheduler.initialize(tree_impl->threadCnt);
		}

		assert(tree_impl->candidates.size() == tree_impl->candidatesMask.size());
		CleanMorphs filterMorphs(tree_impl, tree_impl->candidates.size(), tree_impl->candidates, tree_impl->candidatesMask, newCandidates);
		tbb::parallel_for(
				tbb::blocked_range<size_t>(0, tree_impl->candidates.size()),
				filterMorphs, tbb::auto_partitioner(), tbbCtx);
		tree_impl->candidates.swap(newCandidates);
		tree_impl->candidatesMask.resize(tree_impl->candidates.size());
		tree_impl->candidatesMask.assign(tree_impl->candidates.size(), true);
		assert(tree_impl->candidates.size() == tree_impl->candidatesMask.size());

		newCandidates.clear();
		newCandidates.shrink_to_fit();
	} else {
		throw std::runtime_error("Cannot filter morphs. No tree associated with this instance.");
	}
}



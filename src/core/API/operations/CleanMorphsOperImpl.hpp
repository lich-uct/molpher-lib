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

#ifndef MOLPHER_LIB_CLEANMORPHSOPERIMPL_HPP
#define MOLPHER_LIB_CLEANMORPHSOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class CleanMorphsOper::CleanMorphsOperImpl : public TreeOperation::TreeOperationImpl {
private:

	ConcurrentMolVector newCandidates;

	class CleanMorphs {
	public:

		CleanMorphs(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_pimpl, size_t globalMorphCount,
					 ConcurrentMolVector &morphs, ConcurrentMaskVector &survivors, ConcurrentMolVector &new_candidates);
		void operator()(const tbb::blocked_range<size_t> &r) const;

	private:
		std::shared_ptr<ExplorationTree::ExplorationTreeImpl> mTreePimpl;
		size_t mGlobalMorphCount;
		ConcurrentMolVector &mMorphs;
		ConcurrentMaskVector &mSurvivors;
		ConcurrentMolVector &mNewCandidates;
	};

public:

	CleanMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree);
	CleanMorphsOperImpl();

	void operator()();
};

#endif //MOLPHER_LIB_CLEANMORPHSOPERIMPL_HPP

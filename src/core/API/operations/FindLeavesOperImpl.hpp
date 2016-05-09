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

#ifndef FINDLEAVESOPERIMPL_HPP
#define	FINDLEAVESOPERIMPL_HPP

#include <memory>

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class FindLeavesOper::FindLeavesOperImpl : public TreeOperation::TreeOperationImpl {
    
    private:
        class FindLeaves
        {
            public:
                FindLeaves(ConcurrentMolVector &leaves, bool increment_iters_without_dist_improve);
                void operator()(const TreeMap::range_type &tree) const;

            private:
                ConcurrentMolVector &mLeaves;
                bool mIncrementDistImproveCounter;
        };
    
        ConcurrentMolVector leaves;
        bool incrementDistImproveCounter;
        

    public:
        FindLeavesOperImpl(std::shared_ptr<ExplorationTree> expTree, bool increment_iters_without_dist_improve = false);
        FindLeavesOperImpl(bool increment_iters_without_dist_improve = false);
        void operator()();
        void operator()(ConcurrentMolVector& ret);
        
        std::unique_ptr<MolVector> fetchLeaves();
        void fetchLeaves(ConcurrentMolVector& ret);

};

#endif	/* FINDLEAVESOPERIMPL_HPP */


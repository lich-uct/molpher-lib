/* 
 * File:   FindLeavesOperImpl.hpp
 * Author: sichom
 *
 * Created on January 28, 2016, 4:08 PM
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
        bool mIncrementDistImproveCounter;
        

    public:
        FindLeavesOperImpl(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> expTree, bool increment_iters_without_dist_improve = false);
        FindLeavesOperImpl(bool increment_iters_without_dist_improve = false);
        void operator()();
        
        std::shared_ptr<MolVectorAPI> fetchLeaves();

};

#endif	/* FINDLEAVESOPERIMPL_HPP */


/* 
 * File:   FindLeavesOperImpl.hpp
 * Author: sichom
 *
 * Created on January 28, 2016, 4:08 PM
 */

#ifndef FINDLEAVESOPERIMPL_HPP
#define	FINDLEAVESOPERIMPL_HPP

#include "TreeOperationImpl.hpp"

class FindLeavesOper::FindLeavesOperImpl : public TreeOperation::TreeOperationImpl {
    
//    friend class ExplorationTree;
    
    // TODO fit to new implementation
    
        class FindLeaves
        {
        public:
            FindLeaves(ExplorationTree::MoleculePointerVector &leaves, bool increment_iters_without_dist_improve);
            void operator()(
                const PathFinderContext::CandidateMap::range_type &candidates) const;

        private:
            ExplorationTree::MoleculePointerVector &mLeaves;
            bool mIncrementDistImproveCounter;
        };
    
        ExplorationTree::MoleculePointerVector leaves;
        bool mIncrementDistImproveCounter;
        

    public:
        FindLeavesOper(ExplorationTree& expTree, bool increment_iters_without_dist_improve = false);
        FindLeavesOper();
        void operator()();
        
        const std::vector<MolpherMol>& fetchLeaves();

};

#endif	/* FINDLEAVESOPERIMPL_HPP */


/* 
 * File:   FindLeaves.hpp
 * Author: sichom
 *
 * Created on October 15, 2015, 11:27 AM
 */

#ifndef FINDLEAVES_HPP
#define	FINDLEAVES_HPP

#include "TreeOperation.hpp"

class FindLeavesOper : public TreeOperation {
    
    friend class ExplorationTree;
    
        class FindLeaves
        {
        public:
            FindLeaves(ExplorationTree::MoleculePointerVector &leaves);
            void operator()(
                const PathFinderContext::CandidateMap::range_type &candidates) const;

        private:
            ExplorationTree::MoleculePointerVector &mLeaves;
        };
    
        ExplorationTree::MoleculePointerVector leaves;

    public:
        FindLeavesOper(ExplorationTree& expTree);
        void operator()();

};

#endif	/* FINDLEAVES_HPP */


/* 
 * File:   FindLeaves.hpp
 * Author: sichom
 *
 * Created on October 15, 2015, 11:27 AM
 */

#ifndef FINDLEAVES_HPP
#define	FINDLEAVES_HPP

#include <memory>

#include "TreeOperation.hpp"

class FindLeavesOper : public TreeOperation {
    
//    friend class ExplorationTree;
    
    public:
        class FindLeavesOperImpl;
        
        FindLeavesOper(std::shared_ptr<ExplorationTree> expTree, bool increment_iters_without_dist_improve = false);
        FindLeavesOper(bool increment_iters_without_dist_improve = false);
        void operator()();
        
        std::vector<std::shared_ptr<MolpherMol> > fetchLeaves();
    
    private:
        std::shared_ptr<FindLeavesOperImpl> pimpl;

};

#endif	/* FINDLEAVES_HPP */


/* 
 * File:   PruneTreeOper.hpp
 * Author: sichom
 *
 * Created on October 29, 2015, 2:05 PM
 */

#ifndef PRUNETREEOPER_HPP
#define	PRUNETREEOPER_HPP

#include "TreeOperation.hpp"

class PruneTreeOper : public TreeOperation {
    
    public:
        class PruneTreeOperImpl;
        PruneTreeOper(std::shared_ptr<ExplorationTree> expTree);
        PruneTreeOper();
        virtual void operator()();
        
//        const std::vector<std::string>& getPruned();
        
    private:
        std::shared_ptr<PruneTreeOperImpl> pimpl;
        
};

#endif	/* PRUNETREEOPER_HPP */


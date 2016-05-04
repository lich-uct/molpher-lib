/* 
 * File:   TraverseOper.hpp
 * Author: sichom
 *
 * Created on October 29, 2015, 9:46 AM
 */

#ifndef TRAVERSEOPER_HPP
#define	TRAVERSEOPER_HPP

#include "TreeOperation.hpp"
#include "operations/callbacks/TraverseCallback.hpp"

class TraverseOper : public TreeOperation {
    
    public:
        class TraverseOperImpl;
        
        TraverseOper(std::shared_ptr<ExplorationTree> expTree, TraverseCallback& callback);
        TraverseOper(TraverseCallback& callback);
        TraverseOper(std::shared_ptr<ExplorationTree> expTree, TraverseCallback& callback, const std::string& rootSMILES);
        virtual void operator()();
        
    private:
        std::shared_ptr<TraverseOperImpl> pimpl;
};


#endif	/* TRAVERSEOPER_HPP */


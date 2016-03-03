/* 
 * File:   ExtendTreeOper.hpp
 * Author: sichom
 *
 * Created on October 22, 2015, 5:02 PM
 */

#ifndef EXTENDTREEOPER_HPP
#define	EXTENDTREEOPER_HPP

#include <tbb/parallel_scan.h>

#include "TreeOperation.hpp"

class ExtendTreeOper : public TreeOperation {
    
public:
    class ExtendTreeOperImpl;
    
    ExtendTreeOper(std::shared_ptr<ExplorationTree> expTree);
    ExtendTreeOper();
    virtual void operator()();
    
private:
    std::shared_ptr<ExtendTreeOperImpl> pimpl;
    
};

#endif	/* EXTENDTREEOPER_HPP */


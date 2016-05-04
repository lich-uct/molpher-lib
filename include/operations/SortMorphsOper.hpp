/* 
 * File:   SortMorphsOper.hpp
 * Author: sichom
 *
 * Created on October 21, 2015, 3:14 PM
 */

#ifndef SORTMORPHSOPER_HPP
#define	SORTMORPHSOPER_HPP

#include "TreeOperation.hpp"

class SortMorphsOper : public TreeOperation {
    
public:
    class SortMorphsOperImpl;
    SortMorphsOper(std::shared_ptr<ExplorationTree> expTree);
    SortMorphsOper();
    virtual void operator()();
    
private:
    std::shared_ptr<SortMorphsOperImpl> pimpl;
};

#endif	/* SORTMORPHSOPER_HPP */


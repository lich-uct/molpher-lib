/* 
 * File:   SortMorphsOperImpl.hpp
 * Author: sichom
 *
 * Created on March 1, 2016, 1:16 PM
 */

#ifndef SORTMORPHSOPERIMPL_HPP
#define	SORTMORPHSOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class SortMorphsOper::SortMorphsOperImpl : public TreeOperation::TreeOperationImpl {

private:
    class CompareMorphs {
    public:
        bool operator()(const std::shared_ptr<MolpherMol> &a, const std::shared_ptr<MolpherMol> &b) const;
    };
    
public:
    SortMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree);
    SortMorphsOperImpl();
    virtual void operator()();
};

#endif	/* SORTMORPHSOPERIMPL_HPP */


/* 
 * File:   FilterMorphsOper.hpp
 * Author: sichom
 *
 * Created on October 22, 2015, 10:18 AM
 */

#ifndef FILTERMORPHSOPER_HPP
#define	FILTERMORPHSOPER_HPP

#include "TreeOperation.hpp"

class FilterMorphsOper : public TreeOperation {

public:
    class FilterMorphsOperImpl;

    enum MorphFilters : int {
        PROBABILITY = 1 << 0,
        WEIGHT = 1 << 1,
        SYNTHESIS = 1 << 2,
        MAX_DERIVATIONS = 1 << 3,
        DUPLICATES = 1 << 4,
        HISTORIC_DESCENDENTS = 1 << 5,
        ALL = PROBABILITY | WEIGHT | SYNTHESIS | MAX_DERIVATIONS | DUPLICATES | HISTORIC_DESCENDENTS
    };

    FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, bool verbose = false);
    FilterMorphsOper(bool verbose = false);
    FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, MorphFilters filters, bool verbose = false);
    FilterMorphsOper(MorphFilters filters, bool verbose = false);
    virtual void operator()();
    
private:
    std::shared_ptr<FilterMorphsOperImpl> pimpl;
};


#endif	/* FILTERMORPHSOPER_HPP */


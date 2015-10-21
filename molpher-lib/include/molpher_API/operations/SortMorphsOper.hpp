/* 
 * File:   SortMorphsOper.hpp
 * Author: sichom
 *
 * Created on October 21, 2015, 3:14 PM
 */

#ifndef SORTMORPHSOPER_HPP
#define	SORTMORPHSOPER_HPP

#include "TreeOperation.hpp"

class SortMoprhsOper : public TreeOperation {

private:
    class CompareMorphs {
    public:
        bool operator()(const MolpherMolecule &a, const MolpherMolecule &b) const;
    };
    
public:
    SortMoprhsOper(ExplorationTree& expTree);
    virtual void operator()();
};

#endif	/* SORTMORPHSOPER_HPP */


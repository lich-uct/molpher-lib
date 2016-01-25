/* 
 * File:   TreeOperation.hpp
 * Author: sichom
 *
 * Created on October 13, 2015, 12:43 PM
 */

#ifndef TREEOPERATION_HPP
#define	TREEOPERATION_HPP

#include "../ExplorationTree.hpp"
#include <string>

class ExplorationTree;

class TreeOperation {

protected:
    ExplorationTree* tree;
    int threadCnt;
    
    PathFinderContext& fetchTreeContext();
    ExplorationTree::MoleculeVector& fetchGeneratedMorphs();
    ExplorationTree::BoolVector& fetchGeneratedMorphsMask();
//    void fetchLeaves(ExplorationTree::MoleculePointerVector&);
    
public:
    TreeOperation(ExplorationTree& expTree);
    TreeOperation();
    virtual ~TreeOperation();
    virtual void operator()() = 0;
    
    ExplorationTree* getTree();
    void setTree(ExplorationTree& tree);
};

#endif	/* TREEOPERATION_HPP */


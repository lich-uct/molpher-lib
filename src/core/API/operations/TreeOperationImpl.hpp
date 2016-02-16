/* 
 * File:   TreeOperationImpl.hpp
 * Author: sichom
 *
 * Created on January 27, 2016, 2:06 PM
 */

#ifndef TREEOPERATIONIMPL_HPP
#define	TREEOPERATIONIMPL_HPP

#include <string>

#include "data_structs/ExplorationTree.hpp"
#include "operations/TreeOperation.hpp"

class TreeOperation::TreeOperationImpl {

protected:
    std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_pimpl;
    std::shared_ptr<ExplorationTree> tree;
    
public:
    TreeOperationImpl(std::shared_ptr<ExplorationTree> expTree);
    TreeOperationImpl();
    virtual void operator()();
    
    std::shared_ptr<ExplorationTree> getTree();
    void setTree(std::shared_ptr<ExplorationTree> tree);
};

#endif	/* TREEOPERATIONIMPL_HPP */


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

private:
    std::weak_ptr<ExplorationTree::ExplorationTreeImpl> tree;
    
public:
    TreeOperationImpl(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> expTree);
    TreeOperationImpl();
    virtual void operator()() = 0;
    
    std::shared_ptr<ExplorationTree::ExplorationTreeImpl> getTree();
    void setTree(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree);
};

#endif	/* TREEOPERATIONIMPL_HPP */


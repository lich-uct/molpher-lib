/* 
 * File:   TreeOperation.hpp
 * Author: sichom
 *
 * Created on October 13, 2015, 12:43 PM
 */

#ifndef TREEOPERATION_HPP
#define	TREEOPERATION_HPP

#include <string>
#include <memory>

class ExplorationTree;

class TreeOperation {
    
public:
    class TreeOperationImpl;
    
    TreeOperation(std::shared_ptr<ExplorationTree> expTree);
    TreeOperation();
    virtual ~TreeOperation();
    virtual void operator()() = 0;
    
    std::shared_ptr<ExplorationTree> getTree();
    void setTree(std::shared_ptr<ExplorationTree> tree);
    
protected:
    void setTreeOperPimpl(std::shared_ptr<TreeOperationImpl> impl);
    
private:
    std::shared_ptr<TreeOperationImpl> pimpl;
};

#endif	/* TREEOPERATION_HPP */


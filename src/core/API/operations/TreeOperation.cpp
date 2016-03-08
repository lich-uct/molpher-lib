
#include "operations/TreeOperation.hpp"
#include "TreeOperationImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"


TreeOperation::TreeOperation(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new TreeOperation::TreeOperationImpl(expTree))
{
    // no action
}

TreeOperation::TreeOperation() :
pimpl(new TreeOperation::TreeOperationImpl())
{
    // no action
}

TreeOperation::~TreeOperation() {
    // no action
}

std::shared_ptr<ExplorationTree> TreeOperation::getTree() {
    return pimpl->getTree();
}

void TreeOperation::setTree(std::shared_ptr<ExplorationTree> tree) {
    pimpl->setTree(tree);
}

TreeOperation::TreeOperationImpl::TreeOperationImpl(std::shared_ptr<ExplorationTree> expTree) :
tree(expTree)
{
    // no action
}

TreeOperation::TreeOperationImpl::TreeOperationImpl() {
    // no action
}

std::shared_ptr<ExplorationTree> TreeOperation::TreeOperationImpl::getTree() {
    return tree;
}

void TreeOperation::TreeOperationImpl::setTree(std::shared_ptr<ExplorationTree> tree) {
    this->tree = tree;
}

void TreeOperation::setTreeOperPimpl(std::shared_ptr<TreeOperationImpl> pimpl) {
    this->pimpl = pimpl;
}


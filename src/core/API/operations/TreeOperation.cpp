
#include "operations/TreeOperation.hpp"
#include "TreeOperationImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"


TreeOperation::TreeOperation(ExplorationTree& expTree) : 
pimpl(new TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree.pimpl))
{
    // no action
}

TreeOperation::TreeOperation()
{
    // no action
}

TreeOperation::~TreeOperation() {
    // no action
}

//PathFinderContext& TreeOperation::fetchTreeContext() {
//    return tree->context;
//}
//
//ExplorationTree::MoleculeVector& TreeOperation::fetchGeneratedMorphs() {
//    return tree->candidateMoprhs;
//}
//
//ExplorationTree::BoolVector& TreeOperation::fetchGeneratedMorphsMask() {
//    return tree->candidateMorphsMask;
//}


//void TreeOperation::fetchLeaves(ExplorationTree::MoleculePointerVector& leaves) {
//    tree->fetchLeaves(leaves);
//}

//ExplorationTree* TreeOperation::getTree() {
//    return tree;
//}
//
//void TreeOperation::setTree(ExplorationTree& expTree) {
//    tree = &expTree;
//}

std::shared_ptr<ExplorationTree> TreeOperation::getTree() {
    return tree;
}

void TreeOperation::setTree(std::shared_ptr<ExplorationTree> tree) {
    this->tree = tree;
    pimpl->setTree(tree->pimpl);
}

TreeOperation::TreeOperationImpl::TreeOperationImpl(ExplorationTree::ExplorationTreeImpl& expTree) :
tree(expTree)
{
    // no action
}

TreeOperation::TreeOperationImpl::TreeOperationImpl() {
    // no action
}

std::shared_ptr<ExplorationTree::ExplorationTreeImpl> TreeOperation::TreeOperationImpl::getTree() {
    return tree;
}

void TreeOperation::TreeOperationImpl::setTree(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree) {
    this->tree = tree;
}

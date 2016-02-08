
#include "operations/TreeOperation.hpp"
#include "TreeOperationImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"


TreeOperation::TreeOperation(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new TreeOperation::TreeOperationImpl(expTree->pimpl))
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
    return tree.lock();
}

void TreeOperation::setTree(std::shared_ptr<ExplorationTree> tree) {
    this->tree = tree;
    pimpl->setTree(tree->pimpl);
}

TreeOperation::TreeOperationImpl::TreeOperationImpl(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> expTree) :
tree(expTree)
{
    // no action
}

TreeOperation::TreeOperationImpl::TreeOperationImpl() {
    // no action
}

std::shared_ptr<ExplorationTree::ExplorationTreeImpl> TreeOperation::TreeOperationImpl::getTree() {
    return tree.lock();
}

void TreeOperation::TreeOperationImpl::setTree(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree) {
    this->tree = tree;
}

void TreeOperation::TreeOperationImpl::operator()() {
    // no action
}


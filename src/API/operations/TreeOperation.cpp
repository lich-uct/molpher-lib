
#include "molpher_API/operations/TreeOperation.hpp"


TreeOperation::TreeOperation(ExplorationTree& expTree) : tree(&expTree), threadCnt(expTree.threadCount) {
    // no action
}

TreeOperation::TreeOperation() : tree(nullptr), threadCnt(0) {
    // no action
}

TreeOperation::~TreeOperation() {
    // no action
}

PathFinderContext& TreeOperation::fetchTreeContext() {
    return tree->context;
}

ExplorationTree::MoleculeVector& TreeOperation::fetchGeneratedMorphs() {
    return tree->candidateMoprhs;
}

ExplorationTree::BoolVector& TreeOperation::fetchGeneratedMorphsMask() {
    return tree->candidateMorphsMask;
}


//void TreeOperation::fetchLeaves(ExplorationTree::MoleculePointerVector& leaves) {
//    tree->fetchLeaves(leaves);
//}

ExplorationTree* TreeOperation::getTree() {
    return tree;
}

void TreeOperation::setTree(ExplorationTree& expTree) {
    tree = &expTree;
}


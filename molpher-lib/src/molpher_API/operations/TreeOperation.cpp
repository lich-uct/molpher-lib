
#include "molpher_API/operations/TreeOperation.hpp"


TreeOperation::TreeOperation(ExplorationTree& expTree) : tree(expTree), threadCnt(expTree.threadCount) {
    // no action
}

PathFinderContext& TreeOperation::fetchTreeContext() {
    return tree.context;
}

ExplorationTree::MoleculeVector& TreeOperation::fetchGeneratedMorphs() {
    return tree.candidateMoprhs;
}

ExplorationTree::BoolVector& TreeOperation::fetchGeneratedMorphsMask() {
    return tree.candidateMorphsMask;
}


void TreeOperation::fetchLeaves(ExplorationTree::MoleculeVector& leaves) {
    tree.fetchLeaves(leaves);
}
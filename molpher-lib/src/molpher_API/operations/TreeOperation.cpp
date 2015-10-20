
#include "molpher_API/operations/TreeOperation.hpp"


TreeOperation::TreeOperation(ExplorationTree& expTree) : tree(expTree), threadCnt(expTree.threadCount) {
    // no action
}

PathFinderContext& TreeOperation::fetchTreeContext() {
    return tree.context;
}

ExplorationTree::MoleculeVector& TreeOperation::fetchPutativeLeaves() {
    return tree.putativeLeaves;
}

void TreeOperation::fetchLeaves(ExplorationTree::MoleculeVector& leaves) {
    tree.fetchLeaves(leaves);
}
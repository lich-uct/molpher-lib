
#include "molpher_API/operations/TreeOperation.hpp"


TreeOperation::TreeOperation(ExplorationTree& expTree) : tree(expTree) {
    // no action
}

PathFinderContext& TreeOperation::fetchTreeContext() {
    return tree.context;
}
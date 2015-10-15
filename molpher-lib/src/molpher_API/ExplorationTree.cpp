
#include <stdexcept>

#include "molpher_API/ExplorationTree.hpp"

ExplorationTree::ExplorationTree(IterationSnapshot snp) {
    PathFinderContext::SnapshotToContext(snp, context);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    setParams(params);
}

ExplorationTree::ExplorationTree(const ExplorationParameters& params) {
    setParams(params);
}

void ExplorationTree::setParams(const ExplorationParameters& params) {
    if (params.valid()) {
        IterationSnapshot snp = params.iterSnapshot;
        PathFinderContext::SnapshotToContext(snp, context);
    } else {
        throw std::runtime_error("Exploration parameters are invalid.");
    }
}

ExplorationTreeSnapshot ExplorationTree::createSnapshot() const {
    IterationSnapshot snp;
    PathFinderContext::ContextToSnapshot(context, snp);
    return ExplorationTreeSnapshot(snp);
}

ExplorationTree ExplorationTree::createFromSnapshot(ExplorationTreeSnapshot snapshot) {
    return ExplorationTree(snapshot.iterSnapshot);
}
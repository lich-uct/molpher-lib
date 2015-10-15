
#include <stdexcept>

#include "molpher_API/ExplorationTree.hpp"
#include "molpher_API/operations/FindLeavesOper.hpp"

ExplorationTree::ExplorationTree(IterationSnapshot snp) : threadCount(0) {
    PathFinderContext::SnapshotToContext(snp, context);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : threadCount(0) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    setParams(params);
}

ExplorationTree::ExplorationTree(const ExplorationParameters& params) : threadCount(0) {
    setParams(params);
}

void ExplorationTree::setParams(const ExplorationParameters& params) {
    if (params.valid()) {
        IterationSnapshot snp = params.iterSnapshot;
        PathFinderContext::SnapshotToContext(snp, context);
        PathFinderContext::CandidateMap::accessor ac;
        context.candidates.insert(ac, context.source.smile);
        ac->second = context.source;
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

std::vector<MolpherMol> ExplorationTree::fetchLeaves() {
    FindLeavesOper op(*this);
    op();
    std::vector<MolpherMol> ret;
    for (MoleculeVector::iterator it = op.leaves.begin(); it != op.leaves.end(); it++) {
        ret.push_back(MolpherMol(*it));
    }
    return ret;
}

void ExplorationTree::setThreadCount(int threadCnt) {
    threadCount = threadCnt;
}

int ExplorationTree::getThreadCount() {
    return threadCount;
}
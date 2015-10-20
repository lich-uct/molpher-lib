
#include <stdexcept>
#include <iostream>

#include "molpher_API/ExplorationTree.hpp"
#include "molpher_API/operations/FindLeavesOper.hpp"
#include "molpher_API/operations/PutativeExtendOper.hpp"

ExplorationTree::ExplorationTree(IterationSnapshot& snp) : threadCount(0) {
    PathFinderContext::SnapshotToContext(snp, context);
    initCandidates(snp);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : threadCount(0) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    setParams(params);
}

ExplorationTree::ExplorationTree(ExplorationParameters& params) : threadCount(0) {
    setParams(params);
}

void ExplorationTree::initCandidates(IterationSnapshot& snp) {
    if (snp.target.smile.empty()) {
        snp.target.smile = "C";
        std::cout << "WARNING: No target specified. Inserting default: 'C'" << std::endl;
    }
    PathFinderContext::SnapshotToContext(snp, context);
    PathFinderContext::CandidateMap::accessor ac;
    context.candidates.insert(ac, context.source.smile);
    ac->second = context.source;
}

void ExplorationTree::setParams(ExplorationParameters& params) {
    if (params.valid()) {
        IterationSnapshot snp = params.iterSnapshot;
        initCandidates(snp);
    } else {
        throw std::runtime_error("Exploration parameters are invalid.");
    }
}

ExplorationTreeSnapshot ExplorationTree::createSnapshot() const {
    IterationSnapshot snp;
    PathFinderContext::ContextToSnapshot(context, snp);
    return ExplorationTreeSnapshot(snp);
}

ExplorationTree ExplorationTree::createFromSnapshot(ExplorationTreeSnapshot& snapshot) {
    return ExplorationTree(snapshot.iterSnapshot);
}

void ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret) {
    ExplorationTree::MoleculeVector leaves;
    fetchLeaves(leaves);
    for (MoleculeVector::iterator it = leaves.begin(); it != leaves.end(); it++) {
        ret.push_back(MolpherMol(*it));
    }
}

void ExplorationTree::fetchLeaves(ExplorationTree::MoleculeVector& leaves) {
    FindLeavesOper op(*this);
    op();
    for (auto leaf : op.leaves) {
        leaves.push_back(leaf);
    }
}

void ExplorationTree::putativeExtend() {
    PutativeExtendOper(*this)();
}

void ExplorationTree::setThreadCount(int threadCnt) {
    threadCount = threadCnt;
}

int ExplorationTree::getThreadCount() {
    return threadCount;
}
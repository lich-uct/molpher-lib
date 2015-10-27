
#include <stdexcept>
#include <iostream>

#include "molpher_API/ExplorationTree.hpp"
#include "molpher_API/operations/FindLeavesOper.hpp"
#include "molpher_API/operations/GenerateMorphsOper.hpp"
#include "molpher_API/operations/SortMorphsOper.hpp"
#include "molpher_API/operations/FilterMorphsOper.hpp"
#include "molpher_API/operations/ExtendTreeOper.hpp"

ExplorationTree::ExplorationTree(IterationSnapshot& snp) : threadCount(0) {
    PathFinderContext::SnapshotToContext(snp, context);
    treeInit(snp);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : threadCount(0) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    setParams(params);
}

ExplorationTree::ExplorationTree(ExplorationParameters& params) : threadCount(0) {
    setParams(params);
}

void ExplorationTree::treeInit(IterationSnapshot& snp) {
    if (snp.target.smile.empty()) {
        snp.target.smile = "C";
        std::cerr << "WARNING: No target specified. Inserting default: 'C'" << std::endl;
    }
    PathFinderContext::SnapshotToContext(snp, context);
    PathFinderContext::CandidateMap::accessor ac;
    context.candidates.insert(ac, context.source.smile);
    ac->second = context.source;
}

void ExplorationTree::setParams(ExplorationParameters& params) {
    if (params.valid()) {
        IterationSnapshot snp = params.iterSnapshot;
        treeInit(snp);
    } else {
        throw std::runtime_error("Exploration parameters are invalid.");
    }
}

ExplorationTreeSnapshot* ExplorationTree::createSnapshot() const {
    IterationSnapshot snp;
    PathFinderContext::ContextToSnapshot(context, snp);
    return new ExplorationTreeSnapshot(snp);
}

void ExplorationTree::runOperation(TreeOperation& operation) {
    operation.setTree(*this);
    operation();
}

ExplorationTree* ExplorationTree::createFromSnapshot(ExplorationTreeSnapshot& snapshot) {
    return new ExplorationTree(snapshot.iterSnapshot);
}

void ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret) {
    ExplorationTree::MoleculeVector leaves;
    fetchLeaves(leaves);
    for (MoleculeVector::iterator it = leaves.begin(); it != leaves.end(); it++) {
        ret.push_back(MolpherMol(*it));
    }
}

std::vector<MolpherMol>* ExplorationTree::fetchLeaves() {
    std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
    fetchLeaves(*ret);
    return ret;
}

void ExplorationTree::fetchLeaves(ExplorationTree::MoleculeVector& leaves) {
    FindLeavesOper op(*this);
    op();
    for (auto leaf : op.leaves) {
        leaves.push_back(leaf);
    }
}

void ExplorationTree::generateMorphs() {
    GenerateMoprhsOper(*this)();
}

void ExplorationTree::sortMorphs() {
    SortMoprhsOper(*this)();
}

void ExplorationTree::filterMorphs() {
    FilterMoprhsOper(*this)();
}

void ExplorationTree::filterMorphs(int filters) {
    FilterMoprhsOper(*this, filters)();
}

void ExplorationTree::extend() {
    ExtendTreeOper(*this)();
}


MolpherMol ExplorationTree::fetchMol(const std::string& canonSMILES) {
    PathFinderContext::CandidateMap::accessor ac;
    context.candidates.find(ac, canonSMILES);
    return MolpherMol(ac->second);
}

void ExplorationTree::setThreadCount(int threadCnt) {
    threadCount = threadCnt;
}

int ExplorationTree::getThreadCount() {
    return threadCount;
}

std::vector<MolpherMol> ExplorationTree::getCandidateMorphs() {
    std::vector<MolpherMol> ret;
    for (auto mol : candidateMoprhs) {
        ret.push_back(MolpherMol(mol));
    }
    return ret;
}

std::vector<bool> ExplorationTree::getCandidateMorphsMask() {
    std::vector<bool> ret;
    for (auto status : candidateMorphsMask) {
        ret.push_back(status);
    }
    return ret;
}

void ExplorationTree::setCandidateMorphsMask(std::vector<bool> new_mask) {
    if (new_mask.size() != candidateMorphsMask.size()) {
        throw std::runtime_error("The new mask is not the same length as the old one.");
    }
    for (decltype(candidateMorphsMask.size()) idx = 0; idx != candidateMorphsMask.size(); idx++) {
        candidateMorphsMask[idx] = new_mask[idx];
    }
}



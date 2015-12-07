
#include <stdexcept>
#include <iostream>

#include "molpher_API/ExplorationTree.hpp"
#include "molpher_API/operations/FindLeavesOper.hpp"
#include "molpher_API/operations/GenerateMorphsOper.hpp"
#include "molpher_API/operations/SortMorphsOper.hpp"
#include "molpher_API/operations/FilterMorphsOper.hpp"
#include "molpher_API/operations/ExtendTreeOper.hpp"
#include "molpher_API/operations/PruneTreeOper.hpp"
#include "molpher_API/callbacks/EraseSubtreeCallback.hpp"

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
    ExplorationTree::MoleculePointerVector leaves;
    fetchLeaves(leaves);
    for (auto it: leaves) {
        ret.push_back(MolpherMol(*it));
    }
}

const std::vector<MolpherMol>& ExplorationTree::fetchLeaves() {
    std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
    fetchLeaves(*ret);
    return *ret;
}

void ExplorationTree::fetchLeaves(ExplorationTree::MoleculePointerVector& leaves) {
    // FIXME: this should not increment itersWithoutDistImprovement on molecules in the tree (add option to this method that will toggle that).
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
    FilterMorphsOper(*this)();
}

void ExplorationTree::filterMorphs(int filters) {
    FilterMorphsOper(*this, filters)();
}

void ExplorationTree::extend() {
    ExtendTreeOper(*this)();
}

void ExplorationTree::prune() {
    PruneTreeOper(*this)();
}


MolpherMol* ExplorationTree::fetchMol(const std::string& canonSMILES) {
    if (hasMol(canonSMILES)) {
        PathFinderContext::CandidateMap::accessor ac;
        context.candidates.find(ac, canonSMILES);
        return new MolpherMol(ac->second);
    } else {
        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
    }
}

bool ExplorationTree::hasMol(const std::string& canonSMILES) {
    PathFinderContext::CandidateMap::accessor ac;
    context.candidates.find(ac, canonSMILES);
    return !ac.empty();
}

void ExplorationTree::deleteSubtree(const std::string& canonSMILES) {
    if (hasMol(canonSMILES)) {
        PathFinderContext::CandidateMap::accessor ac;
        context.candidates.find(ac, canonSMILES);
        if (ac->second.parentSmile.empty()) {
            throw std::runtime_error("Deleting root of the tree (" + canonSMILES + ") is not permitted.");
        }
        MolpherMol root(ac->second);
        
        PathFinderContext::CandidateMap::accessor acParent;
        context.candidates.find(acParent, ac->second.parentSmile);
        assert(!acParent.empty());
        acParent->second.descendants.erase(root.getSMILES());
        acParent.release();
        ac.release();
        
        EraseSubtreeCallback cb(context);
        cb.processMorph(root);
    } else {
        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
    }
}

void ExplorationTree::setThreadCount(int threadCnt) {
    threadCount = threadCnt;
}

int ExplorationTree::getThreadCount() {
    return threadCount;
}

const std::vector<MolpherMol>& ExplorationTree::getCandidateMorphs() {
    std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
    for (auto& mol : candidateMoprhs) {
        ret->push_back(MolpherMol(mol));
    }
    return *ret;
}

const std::vector<bool>& ExplorationTree::getCandidateMorphsMask() {
    std::vector<bool>* ret = new std::vector<bool>();
    for (auto status : candidateMorphsMask) {
        ret->push_back(status);
    }
    return *ret;
}

void ExplorationTree::setCandidateMorphsMask(const std::vector<bool>& new_mask) {
    if (new_mask.size() != candidateMorphsMask.size()) {
        throw std::runtime_error("The new mask is not the same length as the old one.");
    }
    for (decltype(candidateMorphsMask.size()) idx = 0; idx != candidateMorphsMask.size(); idx++) {
        candidateMorphsMask[idx] = new_mask[idx];
    }
}

ExplorationParameters& ExplorationTree::getParams() {
    ExplorationParameters* ret = new ExplorationParameters();
    IterationSnapshot snp;
    PathFinderContext::ContextToSnapshot(context, snp);
    ret->iterSnapshot = snp;
    return *ret;
}


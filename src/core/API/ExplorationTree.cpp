
#include <stdexcept>
#include <iostream>

#include "core/misc/iteration_serializer.hpp"

#include "data_structs/ExplorationTree.hpp"
#include "operations/FindLeavesOper.hpp"
#include "operations/GenerateMorphsOper.hpp"
#include "operations/SortMorphsOper.hpp"
#include "operations/FilterMorphsOper.hpp"
#include "operations/ExtendTreeOper.hpp"
#include "operations/PruneTreeOper.hpp"
#include "operations/callbacks/EraseSubtreeCallback.hpp"

ExplorationTree::ExplorationTree(IterationSnapshot& snp) : threadCount(0) {
    treeInit(snp);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : threadCount(0) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    setParams(params);
}

ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) : threadCount(0) {
    ExplorationParameters params;
    params.setSourceMol(sourceMolAsSMILES);
    params.setTargetMol(targetMolAsSMILES);
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
    if (context.candidates.empty()) {
        PathFinderContext::SnapshotToContext(snp, context);
        context.source = molpher::iteration::createMoleculeFromSmile(context.source.smile);
        context.target = molpher::iteration::createMoleculeFromSmile(context.target.smile);
        PathFinderContext::CandidateMap::accessor ac;
        context.candidates.insert(ac, context.source.smile);
        ac->second = context.source;
    } else {
        context.jobId = snp.jobId;
        context.iterIdx = snp.iterIdx;
        context.elapsedSeconds = snp.elapsedSeconds;

        context.fingerprintSelector = (FingerprintSelector) snp.fingerprintSelector;
        context.simCoeffSelector = (SimCoeffSelector) snp.simCoeffSelector;
//        context.dimRedSelector = (DimRedSelector) snp.dimRedSelector;

        context.chemOperSelectors.clear();
        context.chemOperSelectors.resize(snp.chemOperSelectors.size(), (ChemOperSelector) 0);
        for (size_t i = 0; i < snp.chemOperSelectors.size(); ++i) {
            context.chemOperSelectors[i] = (ChemOperSelector) snp.chemOperSelectors[i];
        }

        context.params = snp.params;

        context.source = snp.source;
        context.target = snp.target;
        context.decoys = snp.decoys;
        std::cout << "Parameters changed successfully..." << std::endl;
    }
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

void ExplorationTree::fetchLeaves(ExplorationTree::MoleculePointerVector& leaves, bool increase_dist_improve_counter) {
    FindLeavesOper op(*this, increase_dist_improve_counter);
    op();
    for (auto leaf : op.leaves) {
        leaves.push_back(leaf);
    }
}

void ExplorationTree::generateMorphs() {
    GenerateMorphsOper(*this)();
    filterMorphs(FilterMorphsOper::DUPLICATES);
}

void ExplorationTree::sortMorphs() {
    SortMorphsOper(*this)();
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

int ExplorationTree::getGenerationCount() {
    return context.iterIdx;
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

bool ExplorationTree::isPathFound() {
    PathFinderContext::CandidateMap::const_accessor ac;
    return context.candidates.find(ac, this->context.target.smile);
}



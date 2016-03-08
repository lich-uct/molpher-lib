
#include <stdexcept>
#include <iostream>

//#include "core/misc/iteration_serializer.hpp"
#include "selectors/fingerprint_selectors.h"
#include "selectors/chemoper_selectors.h"
#include "selectors/simcoeff_selectors.h"

#include "data_structs/ExplorationTree.hpp"
#include "ExplorationTreeImpl.h"
#include "core/data_structs/ExplorationDataImpl.hpp"
#include "MolpherMolImpl.hpp"
#include "operations/FindLeavesOper.hpp"
#include "operations/FindLeavesOperImpl.hpp"
#include "core/chem/fingerprintStrategy/MorganFngpr.hpp"
#include "core/misc/inout.h"
#include "core/chem/simCoefStrategy/TanimotoSimCoef.hpp"
#include "operations/SortMorphsOper.hpp"
#include "operations/ExtendTreeOper.hpp"
#include "operations/PruneTreeOper.hpp"
//#include "operations/GenerateMorphsOper.hpp"
//#include "operations/SortMorphsOper.hpp"
//#include "operations/FilterMorphsOper.hpp"
//#include "operations/ExtendTreeOper.hpp"
//#include "operations/PruneTreeOper.hpp"
//#include "operations/callbacks/EraseSubtreeCallback.hpp"

ExplorationTree::ExplorationTree() : pimpl(new ExplorationTree::ExplorationTreeImpl()) {
    // no action
}

//ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : 
//pimpl(new ExplorationTree::ExplorationTreeImpl(sourceMolAsSMILES)) 
//{
//    // no action
//}

//ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) : 
//pimpl(new ExplorationTree::ExplorationTreeImpl(sourceMolAsSMILES, targetMolAsSMILES))
//{
//    // no action
//}

std::shared_ptr<ExplorationTree> ExplorationTree::create(const ExplorationData& data) {
    auto new_tree = std::shared_ptr<ExplorationTree>(new ExplorationTree());
    new_tree->update(data);
    return new_tree;
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) {
    ExplorationData data;
    data.setSource(MolpherMol(sourceMolAsSMILES));
    data.setTarget(MolpherMol(targetMolAsSMILES));
    return create(data);
}

//std::shared_ptr<ExplorationTree> ExplorationTree::create(const std::string& sourceMolAsSMILES) {
//    return std::make_shared<ExplorationTree>(sourceMolAsSMILES);
//}


std::shared_ptr<ExplorationData> ExplorationTree::asData() const {
    return pimpl->asData();
}

void ExplorationTree::update(const ExplorationData& data) {
    pimpl->updateData(data, shared_from_this());
}

bool ExplorationTree::hasMol(const std::string& canonSMILES) {
    return pimpl->hasMol(canonSMILES);
}

bool ExplorationTree::hasMol(std::shared_ptr<MolpherMol> mol) {
    return pimpl->hasMol(mol);
}

std::shared_ptr<MolpherMol> ExplorationTree::fetchMol(const std::string& canonSMILES) {
    return pimpl->fetchMol(canonSMILES);
}

void ExplorationTree::runOperation(TreeOperation& operation) {
    pimpl->runOperation(operation, shared_from_this());
}

MolVector ExplorationTree::fetchLeaves(bool increase_dist_improve_counter) {
    return pimpl->fetchLeaves(shared_from_this(), increase_dist_improve_counter);
}

void ExplorationTree::generateMorphs() {
    pimpl->generateMorphs(shared_from_this());
}

std::vector<std::shared_ptr<MolpherMol> > ExplorationTree::getCandidateMorphs() {
    return pimpl->getCandidateMorphs();
}

std::vector<bool> ExplorationTree::getCandidateMorphsMask() {
    return pimpl->getCandidateMorphsMask();
}

void ExplorationTree::sortMorphs() {
    pimpl->sortMorphs(shared_from_this());
}

void ExplorationTree::filterMorphs(bool verbose_output) {
    pimpl->filterMorphs(shared_from_this(), verbose_output);
}

void ExplorationTree::filterMorphs(FilterMorphsOper::MorphFilters filters, bool verbose_output) {
    pimpl->filterMorphs(filters, shared_from_this(), verbose_output);
}

void ExplorationTree::extend() {
    pimpl->extend(shared_from_this());
}

void ExplorationTree::deleteSubtree(const std::string& canonSMILES, bool descendents_only) {
    pimpl->deleteSubtree(canonSMILES, descendents_only);
}

void ExplorationTree::prune() {
    pimpl->prune(shared_from_this());
}

unsigned ExplorationTree::getGenerationCount() {
    return pimpl->getGenerationCount();
}

bool ExplorationTree::isPathFound() {
    return pimpl->isPathFound();
}

void ExplorationTree::traverse(const std::string& rootSMILES, TraverseCallback& callback) {
    pimpl->traverse(shared_from_this(), rootSMILES, callback);
}

void ExplorationTree::traverse(TraverseCallback& callback) {
    pimpl->traverse(shared_from_this(), callback);
}

// pimpl

ExplorationTree::ExplorationTreeImpl::ExplorationTreeImpl() :
generationCnt(0)
, threadCnt(0)
, fingerprint(FP_MORGAN)
, simCoeff(SC_TANIMOTO)
{
    // no action
}

//ExplorationTree::ExplorationTreeImpl::ExplorationTreeImpl(const std::string& sourceMolAsSMILES) : source(sourceMolAsSMILES)
//{
//    // no action
//}

//ExplorationTree::ExplorationTreeImpl::ExplorationTreeImpl(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) : 
//source(sourceMolAsSMILES)
//, target(targetMolAsSMILES)
//, generationCnt(0)
//, threadCnt(0)
//, fingerprint(FP_MORGAN)
//, simCoeff(SC_TANIMOTO)
//{
//    ExplorationData data;
//    data.setSource(MolpherMol(sourceMolAsSMILES));
//    data.setTarget(MolpherMol(targetMolAsSMILES));
//    
//    try {
//        updateFromData(data);
//    } catch(std::runtime_error& err) {
//        throw std::runtime_error("Failed to create an exploration tree instance from invalid data.");
//    }
//}

//ExplorationTree::ExplorationTreeImpl::ExplorationTreeImpl(ExplorationData &data) {
//    try {
//        updateFromData(data);
//    } catch(std::runtime_error& err) {
//        throw std::runtime_error("Failed to create an exploration tree instance from invalid data.");
//    }
//}

void ExplorationTree::ExplorationTreeImpl::updateData(const ExplorationData& data, std::shared_ptr<ExplorationTree> tree)
{
    if (!data.isValid()) {
        throw std::runtime_error("Supplied exploration data are invalid. " 
                "Cannot update this instance.");
    }
    
    bool is_new_tree = false;
    if (treeMap.empty()) {
        is_new_tree = true;
        auto map = data.getTreeMap();
        for (auto& mol : (*map)) {
            auto added_mol = std::make_shared<MolpherMol>(*(mol.second)); // TODO: a simple move should be more efficient and safe here (http://stackoverflow.com/questions/11711034/stdshared-ptr-of-this)
            added_mol->setOwner(tree);
            treeMap.insert(
                std::make_pair(
                    added_mol->getSMILES()
                    , added_mol
                    )
            );
        }
    } else {
        Cerr("This tree is not empty. Only the morphing parameters will be changed.");
    }
    
    if (is_new_tree) {
        candidates.clear();
        auto cndts = data.getCandidates();
        for (auto& mol : (*cndts)) {
            candidates.push_back(std::make_shared<MolpherMol>(*mol)); // TODO: a simple move should be more efficient and safe here (http://stackoverflow.com/questions/11711034/stdshared-ptr-of-this)
        }
        candidatesMask = data.getCandidatesMask();
                
        morphDerivations.clear();
        for (auto& mol_data : data.getDerivationMap()) {
            morphDerivations.insert(mol_data);
        }
        
        source = MolpherMol(*(data.getSource()));
        
        generationCnt = data.getGenerationCount();
    }
        
    threadCnt = data.getThreadCount();
    chemOpers = data.getChemicalOperators();
    fingerprint = data.getFingerprint();
    simCoeff = data.getSimilarityCoefficient();

    target = MolpherMol(*(data.getTarget()));

    params.cntCandidatesToKeep = data.getCntCandidatesToKeep();
    params.cntCandidatesToKeepMax = data.getCntCandidatesToKeepMax();
    params.itThreshold = data.getItThreshold();
    params.cntMaxMorphs = data.getCntMaxMorphs();
    params.distToTargetDepthSwitch = data.getDistToTargetDepthSwitch();
    params.cntMorphsInDepth = data.getCntMorphsInDepth();
    params.cntMorphs = data.getCntMorphs();
    params.minAcceptableMolecularWeight = data.getMinAcceptableMolecularWeight();
    params.maxAcceptableMolecularWeight = data.getMaxAcceptableMolecularWeight();
    
    auto new_data = this->asData();
    if (!new_data->isValid()) {
        throw std::runtime_error("The tree was created with serious "
                "inconsistencies. Check the 'error_snapshot.xml' file for more details...");
        new_data->save("error_snapshot.xml");
    }
}

std::shared_ptr<ExplorationData> ExplorationTree::ExplorationTreeImpl::asData() const {
    auto data = std::make_shared<ExplorationData>();
    
    for (auto mol : candidates) {
        data->addCandidate(*mol);
    }
    
    data->setCandidatesMask(candidatesMask);
    data->setChemicalOperators(chemOpers);
    data->setFingerprint(fingerprint);
    data->setGenerationCount(generationCnt);
    
    for (auto& item : morphDerivations) {
        data->addToDerivationMap(item.first, item.second);
    }
    
    data->setCntCandidatesToKeep(params.cntCandidatesToKeep);
    data->setCntCandidatesToKeepMax(params.cntCandidatesToKeepMax);
    data->setItThreshold(params.itThreshold);
    data->setCntMaxMorphs(params.cntMaxMorphs);
    data->setDistToTargetDepthSwitch(params.distToTargetDepthSwitch);
    data->setCntMorphsInDepth(params.cntMorphsInDepth);
    data->setCntMorphs(params.cntMorphs);
    data->setMinAcceptableMolecularWeight(params.minAcceptableMolecularWeight);
    data->setMaxAcceptableMolecularWeight(params.maxAcceptableMolecularWeight);
    
    data->setSimilarityCoefficient(simCoeff);
    data->setSource(source);
    data->setTarget(target);
    data->setThreadCount(threadCnt);
    
    for (auto& item : treeMap) {
        data->addToTreeMap(item.first, *(item.second));
    }
    
    return data;
}

std::shared_ptr<MolpherMol> ExplorationTree::ExplorationTreeImpl::fetchMol(const std::string& canonSMILES) {
    if (hasMol(canonSMILES)) {
        TreeMap::accessor ac;
        treeMap.find(ac, canonSMILES);
        return ac->second;
    } else {
        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
    }
}

bool ExplorationTree::ExplorationTreeImpl::hasMol(const std::string& canonSMILES) {
    return static_cast<bool>(treeMap.count(canonSMILES));
}

bool ExplorationTree::ExplorationTreeImpl::hasMol(std::shared_ptr<MolpherMol> mol) {
    TreeMap::accessor ac;
    treeMap.find(ac, mol->getSMILES());
    if (ac.empty()) {
        return false;
    } else {
        return ac->second == mol;
    }
}

void ExplorationTree::ExplorationTreeImpl::runOperation(TreeOperation& operation, std::shared_ptr<ExplorationTree> tree) {
    operation.setTree(tree);
    operation();
}

MolVector ExplorationTree::ExplorationTreeImpl::fetchLeaves(std::shared_ptr<ExplorationTree> tree, bool increase_dist_improve_counter) {
    FindLeavesOper find(increase_dist_improve_counter);
    runOperation(find, tree);
    return find.fetchLeaves();
}

void ExplorationTree::ExplorationTreeImpl::fetchLeaves(std::shared_ptr<ExplorationTree> tree, bool increase_dist_improve_counter, ConcurrentMolVector& ret) {
    FindLeavesOper::FindLeavesOperImpl find(tree, increase_dist_improve_counter);
    find.fetchLeaves(ret);
}

void ExplorationTree::ExplorationTreeImpl::generateMorphs(std::shared_ptr<ExplorationTree> tree) {
    GenerateMorphsOper generate;
    runOperation(generate, tree);
}

MolVector ExplorationTree::ExplorationTreeImpl::getCandidateMorphs() {
    MolVector ret;
    for (auto morph : candidates) {
        ret.push_back(morph);
    }
    return ret;
}

std::vector<bool> ExplorationTree::ExplorationTreeImpl::getCandidateMorphsMask() {
    return candidatesMask;
}

void ExplorationTree::ExplorationTreeImpl::sortMorphs(std::shared_ptr<ExplorationTree> tree) {
    SortMorphsOper sort;
    runOperation(sort, tree);
}

void ExplorationTree::ExplorationTreeImpl::filterMorphs(std::shared_ptr<ExplorationTree> tree, bool verbose_output) {
    FilterMorphsOper filter(verbose_output);
    runOperation(filter, tree);
}

void ExplorationTree::ExplorationTreeImpl::filterMorphs(FilterMorphsOper::MorphFilters filters, std::shared_ptr<ExplorationTree> tree, bool verbose_output) {
    FilterMorphsOper filter(filters, verbose_output);
    runOperation(filter, tree);
}

void ExplorationTree::ExplorationTreeImpl::extend(std::shared_ptr<ExplorationTree> tree) {
    ExtendTreeOper extend;
    runOperation(extend, tree);
}

void ExplorationTree::ExplorationTreeImpl::deleteSubtree(const std::string& canonSMILES, bool descendents_only) {
    if (hasMol(canonSMILES)) {
        TreeMap::accessor acRoot; // root of the subtree to remove
        treeMap.find(acRoot, canonSMILES);
        assert(!acRoot.empty());
        
        if (!descendents_only) {
            if (acRoot->second->getParentSMILES().empty()) {
                throw std::runtime_error("Deleting the absolute root of the tree (" + canonSMILES + ") is not permitted.");
            }

            TreeMap::accessor acRootParent;
            treeMap.find(acRootParent, acRoot->second->getParentSMILES());
            assert(!acRootParent.empty());
            acRootParent->second->removeFromDescendants(acRoot->second->getSMILES());
            
            acRootParent.release();
            acRoot.release();
            
            erase(canonSMILES);
        } else {
            std::set<std::string>::const_iterator it;
            auto descendants = acRoot->second->getDescendants();
            for (it = descendants.begin();
                    it != descendants.end(); it++) {
                acRoot->second->removeFromDescendants(*it);
                erase(*it);
            }
            acRoot->second->setItersWithoutDistImprovement(0);
        }
    } else {
        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
    }
}

void ExplorationTree::ExplorationTreeImpl::erase(const std::string& canonSMILES) {
    std::deque<std::string> toErase;
    toErase.push_back(canonSMILES);

    while (!toErase.empty()) {
        std::string current = toErase.front();
        toErase.pop_front();

        TreeMap::accessor ac;
        treeMap.find(ac, current);
        assert(!ac.empty());

        std::set<std::string>::const_iterator it;
        auto descendants = ac->second->getDescendants();
        for (it = descendants.begin();
                it != descendants.end(); it++) {
            toErase.push_back(*it);
        }

//        pruned.push_back(current);
        auto erased_mol = ac->second;
        treeMap.erase(ac);
        erased_mol->removeFromTree();
    }
}

void ExplorationTree::ExplorationTreeImpl::prune(std::shared_ptr<ExplorationTree> tree) {
    PruneTreeOper prune;
    runOperation(prune, tree);
}

unsigned ExplorationTree::ExplorationTreeImpl::getGenerationCount() {
    return generationCnt;
}

bool ExplorationTree::ExplorationTreeImpl::isPathFound() {
    return hasMol(target.getSMILES());
}

void ExplorationTree::ExplorationTreeImpl::traverse(std::shared_ptr<ExplorationTree> tree, const std::string& rootSMILES, TraverseCallback& callback) {
    TraverseOper traverse(tree, callback, rootSMILES);
    traverse();
}

void ExplorationTree::ExplorationTreeImpl::traverse(std::shared_ptr<ExplorationTree> tree, TraverseCallback& callback) {
    TraverseOper traverse(tree, callback);
    traverse();
}

//std::shared_ptr<ExplorationTree::ExplorationTreeImpl> ExplorationTree::ExplorationTreeImpl::createFromData(ExplorationData& data) {
//    return std::make_shared<ExplorationTree::ExplorationTreeImpl>(data);
//}

//void ExplorationTree::ExplorationTreeImpl::runOperation(std::shared_ptr<TreeOperation::TreeOperationImpl> operation) {
//    
//}

//std::shared_ptr<MolVector> ExplorationTree::ExplorationTreeImpl::fetchLeaves(bool increase_dist_improve_counter) {
//    FindLeavesOper::FindLeavesOperImpl op(std::make_shared<ExplorationTree::ExplorationTreeImpl>(this), (bool) increase_dist_improve_counter); // TODO: create an empty deleter so that the shared pointer doesnt kill the the object pointed to by this
//    op();
//    return op.fetchLeaves();
//}


//ExplorationTree::ExplorationTree(IterationSnapshot& snp) : threadCount(0) {
//    treeInit(snp);
//}
//
//ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES) : threadCount(0) {
//    ExplorationParameters params;
//    params.setSourceMol(sourceMolAsSMILES);
//    setParams(params);
//}
//
//ExplorationTree::ExplorationTree(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) : threadCount(0) {
//    ExplorationParameters params;
//    params.setSourceMol(sourceMolAsSMILES);
//    params.setTargetMol(targetMolAsSMILES);
//    setParams(params);
//}
//
//ExplorationTree::ExplorationTree(ExplorationParameters& params) : threadCount(0) {
//    setParams(params);
//}
//
//void ExplorationTree::treeInit(IterationSnapshot& snp) {
//    if (snp.target.smile.empty()) {
//        snp.target.smile = "C";
//        std::cerr << "WARNING: No target specified. Inserting default: 'C'" << std::endl;
//    }
//    if (context.candidates.empty()) {
//        PathFinderContext::SnapshotToContext(snp, context);
//        context.source = molpher::iteration::createMoleculeFromSmile(context.source.smile);
//        context.target = molpher::iteration::createMoleculeFromSmile(context.target.smile);
//        PathFinderContext::CandidateMap::accessor ac;
//        context.candidates.insert(ac, context.source.smile);
//        ac->second = context.source;
//    } else {
//        context.jobId = snp.jobId;
//        context.iterIdx = snp.iterIdx;
//        context.elapsedSeconds = snp.elapsedSeconds;
//
//        context.fingerprintSelector = (FingerprintSelector) snp.fingerprintSelector;
//        context.simCoeffSelector = (SimCoeffSelector) snp.simCoeffSelector;
////        context.dimRedSelector = (DimRedSelector) snp.dimRedSelector;
//
//        context.chemOperSelectors.clear();
//        context.chemOperSelectors.resize(snp.chemOperSelectors.size(), (ChemOperSelector) 0);
//        for (size_t i = 0; i < snp.chemOperSelectors.size(); ++i) {
//            context.chemOperSelectors[i] = (ChemOperSelector) snp.chemOperSelectors[i];
//        }
//
//        context.params = snp.params;
//
//        context.source = snp.source;
//        context.target = snp.target;
//        context.decoys = snp.decoys;
//        std::cout << "Parameters changed successfully..." << std::endl;
//    }
//}
//
//void ExplorationTree::setParams(ExplorationParameters& params) {
//    if (params.valid()) {
//        IterationSnapshot snp = params.iterSnapshot;
//        treeInit(snp);
//    } else {
//        throw std::runtime_error("Exploration parameters are invalid.");
//    }
//}
//
//ExplorationTreeSnapshot* ExplorationTree::createSnapshot() const {
//    IterationSnapshot snp;
//    PathFinderContext::ContextToSnapshot(context, snp);
//    return new ExplorationTreeSnapshot(snp);
//}
//
//void ExplorationTree::runOperation(TreeOperation& operation) {
//    operation.setTree(*this);
//    operation();
//}
//
//ExplorationTree* ExplorationTree::createFromSnapshot(ExplorationTreeSnapshot& snapshot) {
//    return new ExplorationTree(snapshot.iterSnapshot);
//}
//
//void ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret) {
//    ExplorationTree::MoleculePointerVector leaves;
//    fetchLeaves(leaves);
//    for (auto it: leaves) {
//        ret.push_back(MolpherMol(*it));
//    }
//}
//
//const std::vector<MolpherMol>& ExplorationTree::fetchLeaves() {
//    std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
//    fetchLeaves(*ret);
//    return *ret;
//}
//
//void ExplorationTree::fetchLeaves(ExplorationTree::MoleculePointerVector& leaves, bool increase_dist_improve_counter) {
//    FindLeavesOper op(*this, increase_dist_improve_counter);
//    op();
//    for (auto leaf : op.leaves) {
//        leaves.push_back(leaf);
//    }
//}
//
//void ExplorationTree::generateMorphs() {
//    GenerateMorphsOper(*this)();
//    filterMorphs(FilterMorphsOper::DUPLICATES);
//}
//
//void ExplorationTree::sortMorphs() {
//    SortMorphsOper(*this)();
//}
//
//void ExplorationTree::filterMorphs() {
//    FilterMorphsOper(*this)();
//}
//
//void ExplorationTree::filterMorphs(int filters) {
//    FilterMorphsOper(*this, filters)();
//}
//
//void ExplorationTree::extend() {
//    ExtendTreeOper(*this)();
//}
//
//void ExplorationTree::prune() {
//    PruneTreeOper(*this)();
//}
//
//
//MolpherMol* ExplorationTree::fetchMol(const std::string& canonSMILES) {
//    if (hasMol(canonSMILES)) {
//        PathFinderContext::CandidateMap::accessor ac;
//        context.candidates.find(ac, canonSMILES);
//        return new MolpherMol(ac->second);
//    } else {
//        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
//    }
//}
//
//bool ExplorationTree::hasMol(const std::string& canonSMILES) {
//    PathFinderContext::CandidateMap::accessor ac;
//    context.candidates.find(ac, canonSMILES);
//    return !ac.empty();
//}
//
//void ExplorationTree::deleteSubtree(const std::string& canonSMILES) {
//    if (hasMol(canonSMILES)) {
//        PathFinderContext::CandidateMap::accessor ac;
//        context.candidates.find(ac, canonSMILES);
//        if (ac->second.parentSmile.empty()) {
//            throw std::runtime_error("Deleting root of the tree (" + canonSMILES + ") is not permitted.");
//        }
//        MolpherMol root(ac->second);
//        
//        PathFinderContext::CandidateMap::accessor acParent;
//        context.candidates.find(acParent, ac->second.parentSmile);
//        assert(!acParent.empty());
//        acParent->second.descendants.erase(root.getSMILES());
//        acParent.release();
//        ac.release();
//        
//        EraseSubtreeCallback cb(context);
//        cb.processMorph(root);
//    } else {
//        throw std::runtime_error("Molecule (" + canonSMILES + ") is not present in the tree.");
//    }
//}
//
//void ExplorationTree::setThreadCount(int threadCnt) {
//    threadCount = threadCnt;
//}
//
//int ExplorationTree::getThreadCount() {
//    return threadCount;
//}
//
//int ExplorationTree::getGenerationCount() {
//    return context.iterIdx;
//}
//
//const std::vector<MolpherMol>& ExplorationTree::getCandidateMorphs() {
//    std::vector<MolpherMol>* ret = new std::vector<MolpherMol>();
//    for (auto& mol : candidateMoprhs) {
//        ret->push_back(MolpherMol(mol));
//    }
//    return *ret;
//}
//
//const std::vector<bool>& ExplorationTree::getCandidateMorphsMask() {
//    std::vector<bool>* ret = new std::vector<bool>();
//    for (auto status : candidateMorphsMask) {
//        ret->push_back(status);
//    }
//    return *ret;
//}
//
//void ExplorationTree::setCandidateMorphsMask(const std::vector<bool>& new_mask) {
//    if (new_mask.size() != candidateMorphsMask.size()) {
//        throw std::runtime_error("The new mask is not the same length as the old one.");
//    }
//    for (decltype(candidateMorphsMask.size()) idx = 0; idx != candidateMorphsMask.size(); idx++) {
//        candidateMorphsMask[idx] = new_mask[idx];
//    }
//}
//
//ExplorationParameters& ExplorationTree::getParams() {
//    ExplorationParameters* ret = new ExplorationParameters();
//    IterationSnapshot snp;
//    PathFinderContext::ContextToSnapshot(context, snp);
//    ret->iterSnapshot = snp;
//    return *ret;
//}
//
//bool ExplorationTree::isPathFound() {
//    PathFinderContext::CandidateMap::const_accessor ac;
//    return context.candidates.find(ac, this->context.target.smile);
//}



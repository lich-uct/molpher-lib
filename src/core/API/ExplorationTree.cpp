/*
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdexcept>
#include <iostream>

#include <morphing/operators/AddAtom.hpp>
#include <morphing/operators/RemoveAtom.hpp>
#include <morphing/operators/RemoveBond.hpp>
#include <morphing/operators/MutateAtom.hpp>
#include <morphing/operators/AddBond.hpp>
#include <morphing/operators/InterlayAtom.hpp>
#include <morphing/operators/ContractBond.hpp>
#include <morphing/operators/RerouteBond.hpp>

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

ExplorationTree::ExplorationTree() : pimpl(new ExplorationTree::ExplorationTreeImpl()) {
    // no action
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(const ExplorationData& data) {
    auto new_tree = std::shared_ptr<ExplorationTree>(new ExplorationTree());
    new_tree->update(data);
    return new_tree;
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES) {
    ExplorationData data;
    data.setSource(MolpherMol(sourceMolAsSMILES));
    if (!data.getSource()->isValid()) {
		throw std::runtime_error("source molecule is invalid");
    }
    if (targetMolAsSMILES.empty()) {
        throw std::runtime_error("target molecule SMILES is empty");
    }
    data.setTarget(MolpherMol(targetMolAsSMILES));
    if (!data.getTarget()->isValid()) {
		throw std::runtime_error("target molecule is invalid");
    }
    return create(data);
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(std::shared_ptr<MolpherMol> source, std::shared_ptr<MolpherMol> target) {
    ExplorationData data;
    data.setSource(*source);
    if (target) {
        data.setTarget(*target);
    }
    return create(data);
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(const std::string& filename) {
    auto data = ExplorationData::load(filename);
    return create(*data);
}

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

void ExplorationTree::generateMorphs(const std::vector<std::shared_ptr<MorphCollector> >& collectors) {
    pimpl->generateMorphs(shared_from_this(), collectors);
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

void ExplorationTree::save(const std::string& filename) {
    pimpl->save(filename);
}

void ExplorationTree::setCandidateMorphsMask(const std::vector<bool>& new_mask) {
    pimpl->setCandidateMorphsMask(new_mask);
}

void ExplorationTree::setThreadCount(int threadCnt) {
    pimpl->setThreadCount(threadCnt);
}

int ExplorationTree::getThreadCount() {
    return pimpl->getThreadCount();
}

const std::vector<std::shared_ptr<MorphingOperator>> &ExplorationTree::getMorphingOperators() {
    return pimpl->getMorphingOperators();
}

void ExplorationTree::setMorphingOperators(const std::vector<std::shared_ptr<MorphingOperator>> &operators) {
    pimpl->setMorphingOperators(operators);
}

void ExplorationTree::addMorphingOperator(std::shared_ptr<MorphingOperator> operator_) {
    pimpl->addMorphingOperator(operator_);
}

std::shared_ptr<ExplorationTree> ExplorationTree::create(std::shared_ptr<MolpherMol> source) {
    return ExplorationTree::create(source, nullptr);
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
//        Cerr("This tree has already been initialized. Only the morphing parameters will be changed.");
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
        
        source.reset(new MolpherMol(*(data.getSource())));
        
        generationCnt = data.getGenerationCount();
    }
        
    threadCnt = data.getThreadCount();
    fingerprint = data.getFingerprint();
    simCoeff = data.getSimilarityCoefficient();

    AtomLibrary library(AtomLibrary::getDefaultLibrary()); // TODO: user should be allowed to specify their own library
    if (data.getTarget()) {
        target.reset(new MolpherMol(*(data.getTarget())));
        library = AtomLibrary(*target);
    } else {
        target = nullptr;
    }

    chemOperSelectors = data.getChemicalOperators();
    for (auto val : data.getChemicalOperators()) {
        ChemOperSelector selector = (ChemOperSelector) val;
        switch (selector) {
            case OP_ADD_ATOM:
                chemOpers.push_back(std::make_shared<AddAtom>(library));
                break;
            case OP_REMOVE_ATOM:
                chemOpers.push_back(std::make_shared<RemoveAtom>());
                break;
            case OP_ADD_BOND:
                chemOpers.push_back(std::make_shared<AddBond>());
                break;
            case OP_REMOVE_BOND:
                chemOpers.push_back(std::make_shared<RemoveBond>());
                break;
            case OP_MUTATE_ATOM:
                chemOpers.push_back(std::make_shared<MutateAtom>(library));
                break;
            case OP_INTERLAY_ATOM:
                chemOpers.push_back(std::make_shared<InterlayAtom>(library));
                break;
            case OP_BOND_REROUTE:
                chemOpers.push_back(std::make_shared<RerouteBond>());
                break;
            case OP_BOND_CONTRACTION:
                chemOpers.push_back(std::make_shared<ContractBond>());
                break;
            default:
                break;
        }
    }

    params.cntCandidatesToKeep = data.getCntCandidatesToKeep();
    params.cntCandidatesToKeepMax = data.getCntCandidatesToKeepMax();
    params.itThreshold = data.getItThreshold();
    params.cntMaxMorphs = data.getCntMaxMorphs();
    params.distToTargetDepthSwitch = data.getDistToTargetDepthSwitch();
    params.cntMorphsInDepth = data.getCntMorphsInDepth();
    params.cntMorphs = data.getCntMorphs();
    params.minAcceptableMolecularWeight = data.getMinAcceptableMolecularWeight();
    params.maxAcceptableMolecularWeight = data.getMaxAcceptableMolecularWeight();
	params.sascoreMax = data.getSAScoreMax();
    
    auto new_data = this->asData();
    if (!new_data->isValid()) {
        new_data->save("error_snapshot.xml");
        throw std::runtime_error("The tree was created with serious "
                "inconsistencies. Check the 'error_snapshot.xml' file for more details...");
                // TODO: write some code to be able to find out what went wrong
    }
}

std::shared_ptr<ExplorationData> ExplorationTree::ExplorationTreeImpl::asData() const {
    auto data = std::make_shared<ExplorationData>();
    
    for (auto mol : candidates) {
        data->addCandidate(*mol);
    }
    
    data->setCandidatesMask(candidatesMask);
    data->setChemicalOperators(chemOperSelectors);
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
	data->setSAScoreMax(params.sascoreMax);
    
    data->setSimilarityCoefficient(simCoeff);
    data->setSource(*source);
    if (target) {
        data->setTarget(*target);
    }
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
    return treeMap.count(canonSMILES) > 0;
}

bool ExplorationTree::ExplorationTreeImpl::hasMol(std::shared_ptr<MolpherMol> mol) {
    return hasMol(mol->getSMILES());
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

void ExplorationTree::ExplorationTreeImpl::generateMorphs(std::shared_ptr<ExplorationTree> tree, const std::vector<std::shared_ptr<MorphCollector> >& collectors) {
    GenerateMorphsOper generate(collectors);
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
    DefaultSortCallback callback;
    SortMorphsOper sort(callback);
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
    if (target) {
        return hasMol(target);
    } else {
        return false;
    }
}

void ExplorationTree::ExplorationTreeImpl::traverse(std::shared_ptr<ExplorationTree> tree, const std::string& rootSMILES, TraverseCallback& callback) {
    TraverseOper traverse(tree, callback, rootSMILES);
    traverse();
}

void ExplorationTree::ExplorationTreeImpl::traverse(std::shared_ptr<ExplorationTree> tree, TraverseCallback& callback) {
    TraverseOper traverse(tree, callback);
    traverse();
}

void ExplorationTree::ExplorationTreeImpl::save(const std::string& filename) {
    auto data = asData();
    data->save(filename);
}

void ExplorationTree::ExplorationTreeImpl::setCandidateMorphsMask(const std::vector<bool>& new_mask) {
    candidatesMask = new_mask;
}

void ExplorationTree::ExplorationTreeImpl::setThreadCount(int threadCnt) {
    this->threadCnt = threadCnt;
}

int ExplorationTree::ExplorationTreeImpl::getThreadCount() {
    return threadCnt;
}

const std::vector<std::shared_ptr<MorphingOperator>> &ExplorationTree::ExplorationTreeImpl::getMorphingOperators() {
	return chemOpers;
}

void ExplorationTree::ExplorationTreeImpl::setMorphingOperators(
        const std::vector<std::shared_ptr<MorphingOperator>> &operators) {
    chemOpers = operators;
}

void ExplorationTree::ExplorationTreeImpl::addMorphingOperator(std::shared_ptr<MorphingOperator> operator_) {
    chemOpers.push_back(operator_);
}


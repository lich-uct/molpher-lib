//
// Created by sichom on 13.12.21.
//

#include "ReactionOperatorImpl.hpp"

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <core/misc/SynchRand.h>

ReactionOperator::ReactionOperator() :
		MorphingOperator(),
		pimpl(new ReactionOperatorImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

ReactionOperator::ReactionOperator(const std::string& rxnSmarts) :
		MorphingOperator(),
		pimpl(new ReactionOperatorImpl(rxnSmarts))
{
	setMorphingOperatorPimpl(pimpl);
}

void ReactionOperator::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> ReactionOperator::morph() {
	return pimpl->morph();
}

std::string ReactionOperator::getName() const {
	return pimpl->getName();
}

std::shared_ptr<ReactionOperator> ReactionOperator::create(const std::string &rxnSmarts) {
	return std::make_shared<ReactionOperator>(rxnSmarts);
}

// impl

ReactionOperator::ReactionOperatorImpl::ReactionOperatorImpl() :
	MorphingOperatorImpl()
	, original_rdkit(nullptr)
	, reaction(nullptr)
{
	// no action
}

ReactionOperator::ReactionOperatorImpl::ReactionOperatorImpl(const std::string &rxnSMARTS) :
	ReactionOperatorImpl()
{
	this->rxnSMARTS = rxnSMARTS;
	this->reaction.reset(RDKit::RxnSmartsToChemicalReaction(this->rxnSMARTS));
	if (this->reaction) {
		this->reaction->initReactantMatchers();
	} else {
		std:std::runtime_error("Reaction matchers failed to initialize for SMARTS: " + this->rxnSMARTS);
	}
}

void ReactionOperator::ReactionOperatorImpl::setOriginal(std::shared_ptr<MolpherMol> mol) {
	if (mol) {
		original = mol;
		original_rdkit.reset(mol->asRDMol());
	} else {
		throw std::runtime_error("Invalid reference to original.");
	}
}

std::shared_ptr<MolpherMol> ReactionOperator::ReactionOperatorImpl::morph() {
	if (original_rdkit) {
		boost::shared_ptr<RDKit::ROMol> rd_shared(new RDKit::RWMol(*original_rdkit));
		std::vector<boost::shared_ptr<RDKit::ROMol>> reactant{rd_shared};
		std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>>> products = reaction->runReactants(reactant);
		// FIXME: remove products with transformed locked atoms and lock atoms in the resulting molecule
		if (!products.empty()) {
			int random_product_idx = SynchRand::GetRandomNumber(products.size() - 1);
			auto random_product = products[random_product_idx][0].get();
			std::shared_ptr<MolpherMol> ret(new MolpherMol(random_product));
			const std::string& smiles = ret->getSMILES();
//			if (smiles.find(".") != std::string::npos || smiles.find("*") != std::string::npos) {
//				throw std::runtime_error("Invalid pattern in SMILES found:" + smiles);
//			}
			return ret;
		} else {
			return nullptr;
		}
	} else {
		throw std::runtime_error("No starting molecule set. Set the original molecule to morph first.");
	}
}

std::string ReactionOperator::ReactionOperatorImpl::getName() const {
	return "Reaction: " + rxnSMARTS;
}


//
// Created by sichom on 13.12.21.
//

#include "ReactionOperatorImpl.hpp"

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <core/misc/SynchRand.h>

const std::vector<std::string> ReactionOperator::ReactionOperatorImpl::reactions = std::vector<std::string>({
		"[*;H1,H2,H3:1]>>[*:1]-[C:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[N:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[O:2]",
		"[A;H2,H3:1]>>[A:1]=[C:2]",
		"[A;H2,H3:1]>>[A:1]#[C:2]",
		"[A;H2,H3:1]>>[A:1]=[N:2]",
		"[A;H2,H3:1]>>[A:1]=[O:2]",
		"[A;H2,H3:1]>>[A:1]=[S:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[F:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[S:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[Cl:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[Br:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[I:2]",
		"[*;H1,H2,H3:1]>>[*:1]-[O][P:2](=[O])(=[O])[O]",
		"[N:1]1[C:2][C:3]1>>[N:1]1[C:2][C][N]([C])[C][C:3]1",
		"[A:1]1=[A:2][A:3]=[A:4]1>>[A:1]1=[A:2][C:5]=[C:6][A:3]=[A:4]1",
		"[a:1]1[a:2][a:3][a:4][a:5][a:6][a:7][a:8][a:9][a:10]1>>[a:1]1[a:2][a:3][a:4][a:5][a:6]1.[a:7]1[a:8][a:9][a:10][c][c]1",
		"[*:0][*:1]1[*:2][*:3]1[*:4]>>[*:0][*:2][*:4].[*:2][*:3]",
		"[!#6:1]>>[#6:1].[!#6:2]", // think about this
		"[!#7:1]>>[#7:1].[!#7:2]",
		"[!#8:1]>>[#8:1].[!#8:2]",
		"[*:1]-[*:2]>>[*:1]-[C:3]-[*:2]",
		"[*:1]-[*:2]>>[*:1]-[N:3]-[*:2]",
		"[*:1]-[*:2]>>[*:1]-[O:3]-[*:2]",
		"[*:1]-[*:2]>>[*:1]-[C:3](=[O:4])-[*:2]",
		"[*:1]=[*:2]-[*:3]>>[*:1]-[*:2]=[*:3]",
		"[*:1]=[*:2]>>[*:1]-[*:2]",
		"[*:1]-[*:2]>>[*:1]=[*:2]",
		"[*:1][*:2][*:3]>>[*:1]-[*:3].[*:2]",
		"[*:1]-[*:2]>>[*:1].[*:2]",
		"[*:1]=[*:2]>>[*:1]=[O].[*:2]=[O]",
		"[*:1]~[A:2]~[A:3]~[*:4]>>[*:1]-[*:4].[A:2]~[A:3]",
//		"([*:1][A:2].[A:3][*:4])>>([*:1][A:3].[A:2][*:4])", // these take a lot of time
//		"([*:1][*:2].[*:3][*:4])>>[*:1][*:3].[*:2].[*:4]", // these take a lot of time
//		"([*:1].[A:2])>>([*:1]-[A:2])", // this one takes a lot of time?
		"[a:1]1[a:2][a:3][a:4][a:5][a:6]1>>[a:1]1[a:2][a:3][a:4][O:7]1.[a:5]2[a:6]cccc2",
		"[a:1]1[a:2][a:3][a:4][a:5][a:6]1>>[a:1]1[a:2][a:3][a:4][N:7]1(C).[a:5]2[a:6]cccc2",
		"[a:1]1[a:2][a:3][a:4][*:5]1>>[a:1]1[a:2][a:3][a:4][*:5]c1",
		"[*:1][*:2][*:3][*:4]>>[*:1][*:3][*:2][*:4]",
//		"([*:1]~[*:2].[*:3]~[*:4])>>([*:1].[*:2].[*3].[*4])" // not sure if this one is needed
});

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

std::vector<std::shared_ptr<MorphingOperator>> ReactionOperator::getDefaultOperators() {
	return ReactionOperatorImpl::getDefaultOperators();
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
		throw std::runtime_error("Reaction matchers failed to initialize for SMARTS: " + this->rxnSMARTS);
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

std::shared_ptr<MolpherMol> tryProduct(const std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>>>& products, int idx) {
	auto random_product = products[idx][0].get();
	std::shared_ptr<MolpherMol> ret(new MolpherMol(random_product));
	return ret;
}

std::shared_ptr<MolpherMol> ReactionOperator::ReactionOperatorImpl::morph() {
	if (original_rdkit) {
		boost::shared_ptr<RDKit::ROMol> rd_shared(new RDKit::RWMol(*original_rdkit));
		std::vector<boost::shared_ptr<RDKit::ROMol>> reactant{rd_shared};
		std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>>> products = reaction->runReactants(reactant);
		// FIXME: remove products with transformed locked atoms and lock atoms in the resulting molecule
		if (!products.empty()) {
			std::shared_ptr<MolpherMol> ret = nullptr;
			std::set<int> ignore;
			int random_product_idx = SynchRand::GetRandomNumber(products.size() - 1);
			while(!ret || ignore.find(random_product_idx) != ignore.end()) {
				if (products.size() == ignore.size()) {
					throw std::runtime_error("All reaction products became ignored.");
				}
				try {
					ret = tryProduct(products, random_product_idx);
				} catch (const RDKit::MolSanitizeException &exc) {
					ignore.insert(random_product_idx);
					random_product_idx = SynchRand::GetRandomNumber(products.size() - 1);
				}
			}
//			const std::string& smiles = ret->getSMILES();
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

std::vector<std::shared_ptr<MorphingOperator>> ReactionOperator::ReactionOperatorImpl::getDefaultOperators() {
	std::vector<std::shared_ptr<MorphingOperator>> ret;
	ret.reserve(reactions.size());
	for (const std::string& smarts : reactions) {
		ret.push_back(std::make_shared<ReactionOperator>(smarts));
	}
	return ret;
}


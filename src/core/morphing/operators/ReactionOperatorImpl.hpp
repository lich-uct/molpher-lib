//
// Created by sichom on 13.12.21.
//

#ifndef MOLPHER_LIB_REACTIONOPERATORIMPL_HPP
#define MOLPHER_LIB_REACTIONOPERATORIMPL_HPP

#include <GraphMol/ChemReactions/Reaction.h>
#include "morphing/operators/ReactionOperator.hpp"
#include "MorphingOperatorImpl.hpp"
#include "core/misc/global_types.h"

class ReactionOperator::ReactionOperatorImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::unique_ptr<RDKit::ChemicalReaction> reaction;
	std::string rxnSMARTS;
	static const std::vector<std::string> reactions;

public:
	ReactionOperatorImpl();
	ReactionOperatorImpl(const std::string& rxnSMARTS);

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();
	std::string getName() const;

	static std::vector<std::shared_ptr<MorphingOperator>> getDefaultOperators();
};

#endif //MOLPHER_LIB_REACTIONOPERATORIMPL_HPP

//
// Created by sichom on 13.12.21.
//

#ifndef MOLPHER_LIB_REACTIONOPERATOR_HPP
#define MOLPHER_LIB_REACTIONOPERATOR_HPP

#include "morphing/operators/MorphingOperator.hpp"

class ReactionOperator : public MorphingOperator {
public:
	ReactionOperator();
	ReactionOperator(const std::string& rxnSmarts);

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();
	virtual std::string getName() const;

	static std::vector<std::shared_ptr<MorphingOperator>> getDefaultOperators();

private:
	class ReactionOperatorImpl;
	std::shared_ptr<ReactionOperatorImpl> pimpl;
};

#endif //MOLPHER_LIB_REACTIONOPERATOR_HPP

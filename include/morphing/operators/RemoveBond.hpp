//
// Created by sichom on 11/26/17.
//

#ifndef MOLPHER_LIB_REMOVEBOND_HPP
#define MOLPHER_LIB_REMOVEBOND_HPP

#include "morphing/operators/MorphingOperator.hpp"

class RemoveBond : public MorphingOperator {
public:
	RemoveBond();

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();
	virtual std::string getName() const;

	const std::vector<std::pair<unsigned int, unsigned int>>& getOpenBonds();

private:
	class RemoveBondImpl;
	std::shared_ptr<RemoveBondImpl> pimpl;
};

#endif //MOLPHER_LIB_REMOVEBOND_HPP

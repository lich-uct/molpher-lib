//
// Created by sichom on 11/29/17.
//

#ifndef MOLPHER_LIB_BONDCONTRACTION_HPP
#define MOLPHER_LIB_BONDCONTRACTION_HPP

#include "morphing/operators/MorphingOperator.hpp"

class ContractBond : public MorphingOperator {
public:
	ContractBond();

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();
	virtual std::string getName() const;

	const std::vector<std::pair<unsigned int, unsigned int>>& getOpenBonds();

private:
	class ContractBondImpl;
	std::shared_ptr<ContractBondImpl> pimpl;
};

#endif //MOLPHER_LIB_BONDCONTRACTION_HPP

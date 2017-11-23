//
// Created by sichom on 11/22/17.
//

#ifndef MOLPHER_LIB_ADDBOND_HPP
#define MOLPHER_LIB_ADDBOND_HPP

#include "morphing/operators/MorphingOperator.hpp"

class AddBond : public MorphingOperator {
public:
	AddBond();

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();

private:
	class AddBondImpl;
	std::shared_ptr<AddBondImpl> pimpl;
};

#endif //MOLPHER_LIB_ADDBOND_HPP

//
// Created by sichom on 1/8/18.
//

#ifndef MOLPHER_LIB_REROUTEBOND_HPP
#define MOLPHER_LIB_REROUTEBOND_HPP

#include "morphing/operators/MorphingOperator.hpp"

class RerouteBond : public MorphingOperator {
public:
	RerouteBond();

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();

private:
	class RerouteBondImpl;
	std::shared_ptr<RerouteBondImpl> pimpl;
};

#endif //MOLPHER_LIB_REROUTEBOND_HPP

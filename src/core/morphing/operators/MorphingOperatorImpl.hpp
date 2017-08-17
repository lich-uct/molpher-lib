//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP
#define MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP

#include "morphing/operators/MorphingOperator.hpp"

class MorphingOperator::MorphingOperatorImpl {
protected:
	std::shared_ptr<MolpherMol> original;

public:
	MorphingOperatorImpl();

	std::shared_ptr<MolpherMol> getOriginal();
};

#endif //MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP

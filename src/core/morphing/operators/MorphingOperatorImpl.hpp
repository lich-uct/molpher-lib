//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP
#define MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP

#include "morphing/operators/MorphingOperator.hpp"
#include <core/misc/global_types.h>

class MorphingOperator::MorphingOperatorImpl {
protected:
	std::shared_ptr<MolpherMol> original;

	void writeOriginalLockInfo(std::shared_ptr<MolpherMol>);
	void writeOriginalLockInfo(std::shared_ptr<MolpherMol>, AtomIdx);

public:
	MorphingOperatorImpl();

	std::shared_ptr<MolpherMol> getOriginal();
	virtual void setOriginal(std::shared_ptr<MolpherMol>);
	virtual bool supportsMultiple();
};

#endif //MOLPHER_LIB_MORPHINGOPERATORIMPL_HPP

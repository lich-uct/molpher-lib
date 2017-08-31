//
// Created by sichom on 7/25/17.
//

#include "MorphingOperatorImpl.hpp"

MorphingOperator::MorphingOperator() :
pimpl(new MorphingOperatorImpl())
{
	// no action
}

MorphingOperator::~MorphingOperator() = default;

std::shared_ptr<MolpherMol> MorphingOperator::getOriginal() {
	return pimpl->getOriginal();
}

void MorphingOperator::setMorphingOperatorPimpl(std::shared_ptr<MorphingOperator::MorphingOperatorImpl> pimpl) {
	this->pimpl = pimpl;
}

MorphingOperator::MorphingOperatorImpl::MorphingOperatorImpl() :
original(nullptr)
{
	// no action
}

std::shared_ptr<MolpherMol> MorphingOperator::MorphingOperatorImpl::getOriginal() {
	return original;
}

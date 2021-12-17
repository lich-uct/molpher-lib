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

void MorphingOperator::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

MorphingOperator::MorphingOperatorImpl::MorphingOperatorImpl() :
original(nullptr)
{
	// no action
}

std::shared_ptr<MolpherMol> MorphingOperator::MorphingOperatorImpl::getOriginal() {
	return original;
}

void MorphingOperator::MorphingOperatorImpl::writeOriginalLockInfo(std::shared_ptr<MolpherMol> mol) {
	for (int orig_idx = 0; orig_idx != original->getAtomCount(); orig_idx++) {
		auto orig_atom = original->getAtom(orig_idx);
		if (orig_atom->isLocked()) {
			mol->getAtom(orig_idx)->setLockingMask(orig_atom->getLockingMask());
		}
	}
}

void MorphingOperator::MorphingOperatorImpl::writeOriginalLockInfo(std::shared_ptr<MolpherMol> mol, AtomIdx removed_idx) {
	bool after_removed = false;
	for (int orig_idx = 0; orig_idx != original->getAtomCount(); orig_idx++) {
		if (removed_idx == orig_idx) {
			after_removed = true;
			continue;
		}

		int new_mol_idx = orig_idx;
		if (after_removed) {
			new_mol_idx = orig_idx - 1;
		}

		auto orig_atom = original->getAtom(orig_idx);
		if (orig_atom->isLocked()) {
			mol->getAtom(new_mol_idx)->setLockingMask(orig_atom->getLockingMask());
		}
	}
}

void MorphingOperator::MorphingOperatorImpl::setOriginal(std::shared_ptr<MolpherMol> mol) {
	original = mol;
}

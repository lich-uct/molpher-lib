//
// Created by sichom on 8/28/17.
//

#ifndef MOLPHER_LIB_REMOVEATOMIMPL_HPP
#define MOLPHER_LIB_REMOVEATOMIMPL_HPP

#include "morphing/operators/RemoveAtom.hpp"
#include "MorphingOperatorImpl.hpp"
#include "core/misc/global_types.h"

class RemoveAtom::RemoveAtomImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::vector<AtomIdx> marked_atoms;

public:
	RemoveAtomImpl();

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

	const std::vector<unsigned int>& getMarkedIndices();
	std::vector<std::shared_ptr<MolpherAtom>> getMarkedAtoms();
};

#endif //MOLPHER_LIB_REMOVEATOMIMPL_HPP

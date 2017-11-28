//
// Created by sichom on 11/27/17.
//

#ifndef MOLPHER_LIB_MUTATEATOMIMPL_HPP
#define MOLPHER_LIB_MUTATEATOMIMPL_HPP

#include "morphing/AtomLibrary.hpp"
#include "morphing/operators/MutateAtom.hpp"
#include "MorphingOperatorImpl.hpp"
#include "core/misc/global_types.h"

class MutateAtom::MutateAtomImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::vector<std::vector<MolpherAtom> > replacements;
	AtomLibrary atom_library;

public:
	MutateAtomImpl();
	MutateAtomImpl(const AtomLibrary& atom_library);

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

};

#endif //MOLPHER_LIB_MUTATEATOMIMPL_HPP

//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_ADDATOMIMPL_HPP
#define MOLPHER_LIB_ADDATOMIMPL_HPP

#include "morphing/AtomLibrary.hpp"
#include "morphing/operators/AddAtom.hpp"
#include "MorphingOperatorImpl.hpp"
#include "core/misc/global_types.h"

class AddAtom::AddAtomImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::vector<AtomIdx> open_atoms;
	AtomLibrary atom_library;

public:
	AddAtomImpl();
	AddAtomImpl(const AtomLibrary& atom_library);

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

	const std::vector<unsigned int>& getOpenIndices();
	std::vector<std::shared_ptr<MolpherAtom>> getOpenAtoms();

};

#endif //MOLPHER_LIB_ADDATOMIMPL_HPP

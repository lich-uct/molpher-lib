//
// Created by sichom on 11/27/17.
//

#ifndef MOLPHER_LIB_MUTATEATOM_HPP
#define MOLPHER_LIB_MUTATEATOM_HPP

#include <morphing/AtomLibrary.hpp>
#include "morphing/operators/MorphingOperator.hpp"

class MutateAtom : public MorphingOperator {
public:
	MutateAtom();
	MutateAtom(const AtomLibrary& atom_library);

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();

private:
	class MutateAtomImpl;
	std::shared_ptr<MutateAtomImpl> pimpl;
};

#endif //MOLPHER_LIB_MUTATEATOM_HPP

//
// Created by sichom on 11/28/17.
//

#ifndef MOLPHER_LIB_INTERLAYATOM_HPP
#define MOLPHER_LIB_INTERLAYATOM_HPP

#include <morphing/AtomLibrary.hpp>
#include "morphing/operators/MorphingOperator.hpp"

class InterlayAtom : public MorphingOperator {
public:
	InterlayAtom();
	InterlayAtom(const AtomLibrary& atom_library);

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();

private:
	class InterlayAtomImpl;
	std::shared_ptr<InterlayAtomImpl> pimpl;
};

#endif //MOLPHER_LIB_INTERLAYATOM_HPP

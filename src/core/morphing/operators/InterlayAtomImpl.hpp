//
// Created by sichom on 11/28/17.
//

#ifndef MOLPHER_LIB_INTERLAYATOMIMPL_HPP
#define MOLPHER_LIB_INTERLAYATOMIMPL_HPP

#include "morphing/AtomLibrary.hpp"
#include "morphing/operators/InterlayAtom.hpp"
#include "MorphingOperatorImpl.hpp"
#include "core/misc/global_types.h"

class InterlayAtom::InterlayAtomImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::map<AtomIdx, std::vector<BondIdx> > interlay_candidates;
	std::vector<std::shared_ptr<MolpherAtom>> atoms;
	AtomLibrary atom_library;

public:
	InterlayAtomImpl();
	InterlayAtomImpl(const AtomLibrary& atom_library);

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

};

#endif //MOLPHER_LIB_INTERLAYATOMIMPL_HPP

//
// Created by sichom on 11/22/17.
//

#ifndef MOLPHER_LIB_ADDBONDIMPL_HPP
#define MOLPHER_LIB_ADDBONDIMPL_HPP

#include "core/misc/global_types.h"
#include "morphing/operators/AddBond.hpp"
#include "MorphingOperatorImpl.hpp"

class AddBond::AddBondImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::vector<std::pair<AtomIdx, AtomIdx>> open_bonds;

public:
	AddBondImpl();

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

	const std::vector<std::pair<AtomIdx, AtomIdx>>& getOpenBonds();

};

#endif //MOLPHER_LIB_ADDBONDIMPL_HPP

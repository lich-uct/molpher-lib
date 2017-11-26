//
// Created by sichom on 11/26/17.
//

#ifndef MOLPHER_LIB_REMOVEBONDIMPL_HPP
#define MOLPHER_LIB_REMOVEBONDIMPL_HPP

#include "core/misc/global_types.h"
#include "morphing/operators/RemoveBond.hpp"
#include "MorphingOperatorImpl.hpp"

class RemoveBond::RemoveBondImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::vector<std::pair<AtomIdx, AtomIdx>> open_bonds;
	std::vector<BondIdx> open_bonds_rd;

public:
	RemoveBondImpl();

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

	const std::vector<std::pair<AtomIdx, AtomIdx>>& getOpenBonds();

};

#endif //MOLPHER_LIB_REMOVEBONDIMPL_HPP

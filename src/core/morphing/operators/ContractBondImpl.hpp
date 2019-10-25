//
// Created by sichom on 11/29/17.
//

#ifndef MOLPHER_LIB_BONDCONTRACTIONIMPL_HPP
#define MOLPHER_LIB_BONDCONTRACTIONIMPL_HPP

#include "core/misc/global_types.h"
#include "morphing/operators/ContractBond.hpp"
#include "MorphingOperatorImpl.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

class ContractBond::ContractBondImpl : public MorphingOperator::MorphingOperatorImpl {
private:
	std::vector<std::pair<AtomIdx, AtomIdx>> open_bonds; // first index is the atom that stays
	std::vector<BondIdx> open_bonds_rd;

	void tryToOpenBond(RDKit::Bond *bond, bool swap = false);

public:
	ContractBondImpl();

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

	const std::vector<std::pair<AtomIdx, AtomIdx>>& getOpenBonds();

};

#endif //MOLPHER_LIB_BONDCONTRACTIONIMPL_HPP

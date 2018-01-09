//
// Created by sichom on 1/8/18.
//

#ifndef MOLPHER_LIB_REROUTEBONDIMPL_HPP
#define MOLPHER_LIB_REROUTEBONDIMPL_HPP

#include "core/misc/global_types.h"
#include "morphing/operators/RerouteBond.hpp"
#include "MorphingOperatorImpl.hpp"

class RerouteBond::RerouteBondImpl : public MorphingOperator::MorphingOperatorImpl {
public:
	typedef struct {
		BondIdx bondIdx;
		std::vector<AtomIdx> candidates[2];
	} Candidates;

private:
	std::unique_ptr<RDKit::RWMol> original_rdkit;
	std::vector<Candidates> candidates;
public:
	RerouteBondImpl();

	void setOriginal(std::shared_ptr<MolpherMol> mol);
	std::shared_ptr<MolpherMol> morph();

};

#endif //MOLPHER_LIB_REROUTEBONDIMPL_HPP

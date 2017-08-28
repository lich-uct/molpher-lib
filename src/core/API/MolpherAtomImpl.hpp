//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_MOLPHERATOMIMPL_HPP
#define MOLPHER_LIB_MOLPHERATOMIMPL_HPP

#include <map>
#include "data_structs/MolpherAtom.hpp"

class MolpherAtom::MolpherAtomImpl {

	friend class MolpherAtom; // TODO: implementation should be completely separated -> change architecture

private:
	unsigned int atomic_num;
	std::string symbol;
	int formal_charge;
	double mass;
	int locking_mask;

	static const std::map<MolpherAtom::LockingMask, std::string> locks_map;

public:
	MolpherAtomImpl(RDKit::Atom* rd_atom);

};

#endif //MOLPHER_LIB_MOLPHERATOMIMPL_HPP

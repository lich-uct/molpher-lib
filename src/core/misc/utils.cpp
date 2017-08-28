//
// Created by sichom on 8/28/17.
//

#include "data_structs/MolpherAtom.hpp"
#include "utils.hpp"
#include "global_types.h"
#include "inout.h"

std::map<int, MolpherAtom::LockingMask> parse_atom_locks(const RDKit::ROMol& mol) {
	std::map<int, MolpherAtom::LockingMask> ret;

	RDKit::STR_VECT prop_names = mol.getPropList();
	for (MolpherAtom::LockingMask lock : MolpherAtom::atom_locks) {
		std::string lock_prop_name("MOLPHER_" + MolpherAtom::lockToString(lock));
		if (std::find(prop_names.begin(), prop_names.end(), lock_prop_name) != prop_names.end()) {
			std::vector<std::string> indices;
			split(mol.getProp<std::string>(lock_prop_name), ',', std::back_inserter(indices));
			for (std::string& str_idx : indices) {
				int int_idx = std::stoi(str_idx) - 1;
				if (ret.find(int_idx) == ret.end()) {
					ret.insert(std::pair<int, MolpherAtom::LockingMask>(int_idx, lock));
				} else {
					ret.find(int_idx)->second = (MolpherAtom::LockingMask) (ret.find(int_idx)->second | lock);
				}
			}
			mol.clearProp(lock_prop_name);
		}
	}

	return ret;
};
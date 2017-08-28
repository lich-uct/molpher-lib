//
// Created by sichom on 8/28/17.
//

#ifndef MOLPHER_LIB_UTILS_HPP
#define MOLPHER_LIB_UTILS_HPP

#include <map>
#include <string>
#include <vector>
#include <GraphMol/ROMol.h>

std::map<int, MolpherAtom::LockingMask> parse_atom_locks(const RDKit::ROMol& mol);

#endif //MOLPHER_LIB_UTILS_HPP

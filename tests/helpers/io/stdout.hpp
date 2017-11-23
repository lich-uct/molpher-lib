//
// Created by sichom on 6/6/17.
//

#ifndef MOLPHER_LIB_STDOUT_HPP
#define MOLPHER_LIB_STDOUT_HPP

#include <GraphMol/ROMol.h>
#include <data_structs/MolpherMol.hpp>

void print_mol_info(RDKit::ROMol* mol);
void print_morphs(const MolpherMol& parent, const std::vector<std::shared_ptr<MolpherMol>>& mols);
void print_locks(const MolpherMol& mol);

void print(const std::string& text);

#endif //MOLPHER_LIB_STDOUT_HPP

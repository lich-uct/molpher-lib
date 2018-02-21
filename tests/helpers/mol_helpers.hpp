//
// Created by sichom on 6/9/17.
//

#ifndef MOLPHER_LIB_MOL_HELPERS_HPP
#define MOLPHER_LIB_MOL_HELPERS_HPP

#include <string>
#include <data_structs/MolpherMol.hpp>

bool match_substr(std::string smiles_mol, std::string smiles_query);
void print_lock_info(std::shared_ptr<MolpherMol> mol);

#endif //MOLPHER_LIB_MOL_HELPERS_HPP

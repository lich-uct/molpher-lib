/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Vladimir Fiklik
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <string>
#include <vector>
#include <map>

#include <GraphMol/GraphMol.h>

#include "global_types.h"
//#include "data_structs/IterationSnapshot.h"
//#include "data_structs/MolpherMolecule.h"

void SynchCout(const std::string &s);
void SynchCerr(const std::string &s, const std::string prefix = "WARNING: ");

void Cout(const std::string &s);
void Cerr(const std::string &s, const std::string prefix = "WARNING: ");

template<typename Number>
std::string parseNumber(Number num) {
    std::stringstream ss;
    ss << num;
    return ss.str();
}

void WriteRWMolsToSDF(const std::string &file,
    std::vector<RDKit::RWMol *> &mols);
void ReadRWMolsFromSDF(const std::string &file,
    std::vector<RDKit::RWMol *> &mols); // caller is responsible for deallocation

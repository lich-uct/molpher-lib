/*
 Copyright (c) 2012 Peter Szepe
 Copyright (c) 2012 Petr Koupy

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

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/Atom.h>

#include "core/data_structs/MolpherAtom.h"
#include "core/data_structs/MolpherMolData.hpp"
#include "data_structs/MolpherMol.hpp"

//typedef unsigned int JobId;
//typedef unsigned int IterIdx;

typedef ExplicitBitVect Fingerprint;
typedef unsigned int AtomIdx;
typedef unsigned int BondIdx;

typedef unsigned int MolpherAtomIdx;

typedef tbb::concurrent_hash_map<std::string, std::shared_ptr<MolpherMol>> TreeMap;
typedef tbb::concurrent_hash_map<std::string, unsigned int> MorphDerivationMap;
typedef tbb::concurrent_vector<std::shared_ptr<MolpherMol>> ConcurrentMolVector;
typedef std::vector<std::shared_ptr<MolpherMol>> MolVector;

typedef std::vector<MolpherMolData> CandidatesVectorData;
typedef std::vector<bool> CandidatesMaskVectorData;
typedef std::map<std::string, MolpherMolData> TreeMapData;
typedef std::map<std::string, unsigned> MorphDerivationMapData;

template<typename Content>
void concurrent_vector_to_vector(const tbb::concurrent_vector<Content>& in, std::vector<Content>& out) {
    for (auto& item : in) {
        out.push_back(item);
    }
}

template<typename Content>
void vector_to_concurrent_vector(const std::vector<Content>& in, tbb::concurrent_vector<Content>& out) {
    for (auto& item : in) {
        out.push_back(item);
    }
}

/*
 Copyright (c) 2012 Peter Szepe
 Copyright (c) 2012 Petr Koupy
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

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Atom.h>

#include "core/data_structs/MolpherMolData.hpp"
#include "data_structs/MolpherMol.hpp"

typedef ExplicitBitVect Fingerprint;
typedef unsigned int AtomIdx;
typedef unsigned int BondIdx;
typedef int AtomicNum;

typedef tbb::concurrent_hash_map<std::string, std::shared_ptr<MolpherMol>> TreeMap;
typedef tbb::concurrent_hash_map<std::string, unsigned int> MorphDerivationMap;
typedef tbb::concurrent_hash_map<std::string, bool /*dummy*/> ConcurrentSmileSet;
typedef tbb::concurrent_vector<std::string> ConcurrentSmileVector;
typedef tbb::concurrent_vector<bool> ConcurrentMaskVector;
typedef tbb::concurrent_vector<std::shared_ptr<MolpherMol>> ConcurrentMolVector;
typedef std::vector<std::shared_ptr<MolpherMol>> MolVector;

typedef std::vector<MolpherMolData> CandidatesVectorData;
typedef std::vector<bool> CandidatesMaskVectorData;
typedef std::map<std::string, MolpherMolData> TreeMapData;
typedef std::map<std::string, unsigned> MorphDerivationMapData;

template<typename Content>
void concurrent_vector_to_vector(const tbb::concurrent_vector<Content>& in, std::vector<Content>& out) {
    out.clear();
	out.reserve(in.size());
	for (Content item : in) {
        out.push_back(item);
    }
}

template<typename Content>
void vector_to_concurrent_vector(const std::vector<Content>& in, tbb::concurrent_vector<Content>& out) {
	out.clear();
	out.reserve(in.size());
	for (Content item : in) {
        out.push_back(item);
    }
}

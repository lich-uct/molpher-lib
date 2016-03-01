/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Vladimir Fiklik

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
void SynchCerr(const std::string &s, const std::string prefix = "WARNIG: ");

void Cout(const std::string &s);
void Cerr(const std::string &s, const std::string prefix = "WARNIG: ");

template<typename Number>
std::string parseNumber(Number num) {
    std::stringstream ss;
    ss << num;
    return ss.str();
}

//IterationSnapshot Materialize(IterationSnapshotProxy &proxy);

//std::string GenerateFilename(std::string &base,
//    JobId jobId, unsigned int iterIdx);
//std::string GenerateFilename(std::string &base,
//    JobId jobId, std::string name);
//std::string GenerateDirname(std::string &base, JobId jobId);
//
//void WriteMolpherPath(const std::string &file, const std::string &targetSmile,
//    const IterationSnapshot::CandidateMap &candidates);
//
//void WriteSnapshotToFile(const std::string &file, const IterationSnapshot &snp);
//bool ReadSnapshotFromFile(const std::string &file, IterationSnapshot &snp);
//
//void GatherMolphMols(const IterationSnapshot::CandidateMap &toGather,
//    std::map<std::string, MolpherMolecule> &gathered);
//void GatherMolphMols(const std::vector<MolpherMolecule> &toGather,
//    std::map<std::string, MolpherMolecule> &gathered);
//
//void WriteMolphMolsToSDF(const std::string &file,
//    const std::map<std::string, MolpherMolecule> &mols);
//void WriteMolphMolsToSDF(const std::string &file,
//    const std::vector<MolpherMolecule> &mols);
//void ReadMolphMolsFromFile(const std::string &file,
//    std::vector<MolpherMolecule> &mols);

void WriteRWMolsToSDF(const std::string &file,
    std::vector<RDKit::RWMol *> &mols);
void ReadRWMolsFromSDF(const std::string &file,
    std::vector<RDKit::RWMol *> &mols); // caller is responsible for deallocation

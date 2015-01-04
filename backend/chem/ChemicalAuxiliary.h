/*
 Copyright (c) 2012 Peter Szepe

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

#include <vector>

#include <GraphMol/GraphMol.h>

#include "global_types.h"

void SetFormalCharge(int charge, RDKit::Atom &atom);

void SetFormalCharge(int charge, RDKit::RWMol &mol);

void GetAtomTypesFromMol(RDKit::ROMol &mol, std::vector<MolpherAtom> &atoms);

MolpherAtomIdx GetRandomAtom(const std::vector<MolpherAtom> &atoms, RDKit::Atom &atom);

RDKit::Bond *GetRandomNonSingleBond(RDKit::Atom &atom);

bool HasNonSingleBond(RDKit::Atom &atom);

void SetBondOrder(RDKit::Bond &bond, int bondOrder);

void DecreaseBondOrder(RDKit::Bond &bond);

void IncreaseBondOrder(RDKit::Bond &bond);

unsigned int CntFreeOxygens(RDKit::Atom &atom);

int GetMaxBondsMod(AtomicNum atomicNum);

int GetMaxBondsMod(MolpherAtom &atom);

int GetMaxBondsMod(RDKit::Atom &atom);

int CntFreeBonds(RDKit::Atom &atom);

void GetPossibleBondingAtoms(
    RDKit::ROMol &mol, int bondOrder, std::vector<AtomIdx> &bondingAtoms);

void GetAtomsNotInRing(RDKit::ROMol &mol, std::vector<RDKit::Atom *> &atomsNIR);

void GetAtomsWithNotMaxValence(RDKit::ROMol &mol, std::vector<RDKit::Atom *> &atomsNMV);

void CopyMol(RDKit::ROMol &mol, RDKit::RWMol &copy);

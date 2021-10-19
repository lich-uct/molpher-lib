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

#include <map>
#include <queue>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

#include "core/misc/SynchRand.h"
#include "ChemicalAuxiliary.h"

RDKit::Atom* GetRandomAtom(const AtomLibrary& atoms)
{
    const MolpherAtom& random_atm = atoms.getRandomAtom();
	auto* atom = new RDKit::Atom();
    atom->setAtomicNum((int) random_atm.getAtomicNum());
    atom->setFormalCharge((int) random_atm.getFormalCharge());
  	return atom;
}

RDKit::Bond *GetRandomNonSingleBond(RDKit::Atom &atom)
{
    std::vector<RDKit::Bond *> candidates;

    RDKit::Bond *bond;
    RDKit::ROMol &mol = atom.getOwningMol();
    RDKit::ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(&atom);
    while (beg != end) {
        bond = mol[*beg++];
        int bo = RDKit::queryBondOrder(bond);
        if (bo > 1 && bo <= 6) {
            candidates.push_back(bond);
        }
    }

    return candidates[SynchRand::GetRandomNumber(0, candidates.size() - 1)];
}

bool HasNonSingleBond(RDKit::Atom &atom)
{
    RDKit::Bond *bond;
    RDKit::ROMol &mol = atom.getOwningMol();
    RDKit::ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(&atom);
    while (beg != end) {
        bond = mol[*beg++];
        int bo = RDKit::queryBondOrder(bond);
        if (bo > 1 && bo <= 6) {
            return true;
        }
    }
    return false;
}

void DecreaseBondOrder(RDKit::Bond &bond)
{
    int bo = RDKit::queryBondOrder(&bond);
    int newBo;

    if (bo >= 2 && bo <= 6) {
        newBo = bo - 1;
    } else if (bo >= 7 && bo <= 11) {
        newBo = bo - 6;
    } else {
        newBo = 1;
    }
    bond.setBondType(static_cast<RDKit::Bond::BondType>(newBo));
}

void IncreaseBondOrder(RDKit::Bond &bond)
{
    int bo = RDKit::queryBondOrder(&bond);
    int newBo;

    if (bo >= 1 && bo <= 5) {
        newBo = bo + 1;
    } else if (bo >= 7 && bo <= 11) {
        newBo = bo - 5;
    } else {
        newBo = 2;
    }
    bond.setBondType(static_cast<RDKit::Bond::BondType>(newBo));
}

int GetMaxBondsMod(unsigned int atomicNum)
{
    return RDKit::PeriodicTable::getTable()->getDefaultValence(atomicNum);
}

int GetMaxBondsMod(MolpherAtom &atom)
{
    int defaultValence =
        RDKit::PeriodicTable::getTable()->getDefaultValence(atom.getAtomicNum());

    int result = defaultValence + atom.getFormalCharge();

    return result;
}

int GetMaxBondsMod(RDKit::Atom &atom)
{
    int defaultValence =
        RDKit::PeriodicTable::getTable()->getDefaultValence(atom.getAtomicNum());

    int result = defaultValence + atom.getFormalCharge();

    return result;
}

int CntFreeBonds(RDKit::Atom &atom)
{
    int freeBonds = atom.getImplicitValence();
    return freeBonds;
}

void GetPossibleBondingAtoms(
        RDKit::ROMol &mol, int bondOrder, const std::set<int>& fixed_ids, std::vector<AtomIdx> &bondingAtoms)
{
    RDKit::Atom *atom;
    RDKit::ROMol::AtomIterator iter;
    for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
        atom = *iter;
        int cntFreeBonds = CntFreeBonds(*atom);

        if (cntFreeBonds >= bondOrder && fixed_ids.find(atom->getIdx()) == fixed_ids.end()) {
            bondingAtoms.push_back(atom->getIdx());
        }
    }
}

void GetAtomsNotInRing(RDKit::ROMol &mol, std::vector<RDKit::Atom *> &atomsNIR)
{
    RDKit::Atom *atom;
    RDKit::ROMol::AtomIterator iter;
    for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
        atom = *iter;
        if (!RDKit::queryIsAtomInRing(atom)) {
            atomsNIR.push_back(atom);
        }
    }
}

void GetAtomsWithNotMaxValence(RDKit::ROMol &mol, std::vector<RDKit::Atom *> &atomsNMV)
{
    RDKit::Atom *atom;
    RDKit::ROMol::AtomIterator iter;
    for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
        atom = *iter;
        if (CntFreeBonds(*atom)) {
            atomsNMV.push_back(atom);
        }
    }
}

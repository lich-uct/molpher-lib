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

#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"

#include "OpInterlayAtom.hpp"

void OpInterlayAtom::Morph(MorphingData &data, RDKit::RWMol **nMol)
{
    *nMol = new RDKit::RWMol(data.mol);
    RDKit::RWMol *newMol = *nMol;

    RDKit::Atom atom;
    AtomIdx idx = GetRandomAtom(data.atoms, atom);

    if (data.interlayAtomCandidates.find(idx) == data.interlayAtomCandidates.end()) {
        delete newMol;
        *nMol = NULL;
        return;
    }
    if (data.interlayAtomCandidates[idx].size() == 0) {
        delete newMol;
        *nMol = NULL;
        return;
    }

    int randPos =
        SynchRand::GetRandomNumber(data.interlayAtomCandidates[idx].size() - 1);
    RDKit::Bond *bond =
        newMol->getBondWithIdx((data.interlayAtomCandidates[idx])[randPos]);

    AtomIdx beg = bond->getBeginAtomIdx();
    AtomIdx end = bond->getEndAtomIdx();
    RDKit::Bond::BondType bt = bond->getBondType();

    AtomIdx newAtomIdx = newMol->addAtom(&atom); // atom is copied

    newMol->removeBond(beg, end);
    newMol->addBond(beg, newAtomIdx, bt);
    newMol->addBond(newAtomIdx, end, bt);
}

ChemOperSelector OpInterlayAtom::GetSelector()
{
    return OP_INTERLAY_ATOM;
}
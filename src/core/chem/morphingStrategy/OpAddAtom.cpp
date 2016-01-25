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

#include "auxiliary/SynchRand.h"
#include "chem/ChemicalAuxiliary.h"

#include "OpAddAtom.hpp"

void OpAddAtom::Morph(MorphingData &data, RDKit::RWMol **nMol)
{
    *nMol = new RDKit::RWMol(data.mol);
    RDKit::RWMol *newMol = *nMol;

    RDKit::Atom atom;
    GetRandomAtom(data.atoms, atom);

    if (data.addAtomCandidates.size() == 0) {
        delete newMol;
        *nMol = NULL;
        return;
    }

    AtomIdx bindingAtomIdx =
        data.addAtomCandidates[SynchRand::GetRandomNumber(
        data.addAtomCandidates.size() - 1)];

    RDKit::Atom *bindingAtom = newMol->getAtomWithIdx(bindingAtomIdx);

    AtomIdx newAtomIdx = newMol->addAtom(&atom); // atom is copied

    if (HasNonSingleBond(*bindingAtom) && (SynchRand::GetRandomNumber(0, 1) > 0)) {
        RDKit::Bond *bond = GetRandomNonSingleBond(*bindingAtom);
        DecreaseBondOrder(*bond);
    }

    newMol->addBond(bindingAtom->getIdx(), newAtomIdx, RDKit::Bond::SINGLE);
}

ChemOperSelector OpAddAtom::GetSelector()
{
    return OP_ADD_ATOM;
}
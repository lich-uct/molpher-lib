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

#include "OpMutateAtom.hpp"

void OpMutateAtom::Morph(MorphingData &data, RDKit::RWMol **nMol)
{
    *nMol = new RDKit::RWMol(data.mol);
    RDKit::RWMol *newMol = *nMol;

    int randPos = SynchRand::GetRandomNumber(newMol->getNumAtoms() - 1);

    if(data.mutateAtomCandidates[randPos].size() == 0) {
        delete newMol;
        *nMol = NULL;
        return;
    }

    RDKit::Atom atom;
    GetRandomAtom(data.mutateAtomCandidates[randPos], atom);

    newMol->replaceAtom(randPos, &atom); // atom is copied
}

ChemOperSelector OpMutateAtom::GetSelector()
{
    return OP_MUTATE_ATOM;
}
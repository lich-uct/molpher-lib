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
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"

#include "OpRemoveAtom.hpp"

void OpRemoveAtom::Morph(MorphingData &data, RDKit::RWMol **nMol)
{
    *nMol = new RDKit::RWMol(data.mol);
    RDKit::RWMol *newMol = *nMol;

    if (data.removeAtomCandidates.size() == 0) {
        delete newMol;
        *nMol = NULL;
        return;
    }

    // remove a random boundary atom
    int randPos = SynchRand::GetRandomNumber(data.removeAtomCandidates.size() - 1);
    AtomIdx atomIdx = data.removeAtomCandidates[randPos];

    newMol->removeAtom(atomIdx);
    tbb::concurrent_hash_map<std::string,int>::accessor ac;
    data.removed_atoms.insert(ac, RDKit::MolToSmiles(*newMol));
    ac->second = atomIdx;
    ac.release();
}

ChemOperSelector OpRemoveAtom::GetSelector()
{
    return OP_REMOVE_ATOM;
}
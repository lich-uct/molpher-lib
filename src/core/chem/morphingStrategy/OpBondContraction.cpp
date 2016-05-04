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

#include "OpBondContraction.hpp"

void OpBondContraction::Morph(MorphingData &data, RDKit::RWMol **nMol)
{
    *nMol = new RDKit::RWMol(data.mol);
    RDKit::RWMol *newMol = *nMol;

    if (data.bondContractionCandidates.size() == 0) {
        delete newMol;
        *nMol = NULL;
        return;
    }

    int randPos =
        SynchRand::GetRandomNumber(data.bondContractionCandidates.size() - 1);
    RDKit::Bond *bond =
        newMol->getBondWithIdx(data.bondContractionCandidates[randPos]);
    RDKit::Atom *atomToRemove = bond->getEndAtom();
    RDKit::Atom *atomToStay = bond->getBeginAtom();

    std::vector<RDKit::Bond *> bondsToRemove;
    RDKit::Bond *bondToChange;
    RDKit::ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = newMol->getAtomBonds(atomToRemove);
    while (beg != end) {
        bondToChange = (**nMol)[*beg++].get();
        if (bondToChange->getIdx() == bond->getIdx()) {
            continue;
        }
        if (bondToChange->getBeginAtomIdx() == atomToRemove->getIdx()) {
            if (!newMol->getBondBetweenAtoms(
                        atomToStay->getIdx(), bondToChange->getEndAtomIdx()) &&
                    !newMol->getBondBetweenAtoms(
                        bondToChange->getEndAtomIdx(), atomToStay->getIdx())) {
                newMol->addBond(atomToStay->getIdx(),
                    bondToChange->getEndAtomIdx(), bondToChange->getBondType());
            }
            bondsToRemove.push_back(bondToChange);
        }
        if (bondToChange->getEndAtomIdx() == atomToRemove->getIdx()) {
            if (!newMol->getBondBetweenAtoms(
                        bondToChange->getBeginAtomIdx(), atomToStay->getIdx()) &&
                    !newMol->getBondBetweenAtoms(
                        atomToStay->getIdx(), bondToChange->getBeginAtomIdx())) {
                newMol->addBond(bondToChange->getBeginAtomIdx(),
                    atomToStay->getIdx(), bondToChange->getBondType());
            }
            bondsToRemove.push_back(bondToChange);
        }
    }

    AtomIdx beginAtomIdx;
    AtomIdx endAtomIdx;
    for (int i = 0; i < bondsToRemove.size(); ++i) {
        beginAtomIdx = bondsToRemove[i]->getBeginAtomIdx();
        endAtomIdx = bondsToRemove[i]->getEndAtomIdx();
        newMol->removeBond(beginAtomIdx, endAtomIdx);
    }

    newMol->removeAtom(atomToRemove);
}

ChemOperSelector OpBondContraction::GetSelector()
{
    return OP_BOND_CONTRACTION;
}

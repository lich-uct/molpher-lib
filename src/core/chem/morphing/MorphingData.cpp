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

#include <queue>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

#include "auxiliary/SynchRand.h"
#include "chemoper_selectors.h"
#include "chem/ChemicalAuxiliary.h"
#include "MorphingData.h"

MorphingData::MorphingData(
    RDKit::ROMol &molecule,
    RDKit::ROMol &target,
    std::vector<ChemOperSelector> &operators
    ) :
    mol(molecule),
    operators(operators)
{
    GetAtomTypesFromMol(target, atoms);

    for (int i = 0; i < operators.size(); ++i) {
        switch (operators[i]) {
        case OP_ADD_ATOM:
            InitAddAtom();
            break;
        case OP_REMOVE_ATOM:
            InitRemoveAtom();
            break;
        case OP_ADD_BOND:
            InitAddBond();
            break;
        case OP_REMOVE_BOND:
            InitRemoveBond();
            break;
        case OP_MUTATE_ATOM:
            InitMutateAtom();
            break;
        case OP_INTERLAY_ATOM:
            InitInterlayAtom();
            break;
        case OP_BOND_REROUTE:
            InitBondReroute();
            break;
        case OP_BOND_CONTRACTION:
            InitBondContraction();
            break;
        }
    }
}

MorphingData::~MorphingData()
{
    // no-op
}

void MorphingData::InitAddAtom()
{
    int bondOrder = 1;
    GetPossibleBondingAtoms(mol, bondOrder, addAtomCandidates);
}

void MorphingData::InitAddBond()
{
    // Calculate atoms product.

    std::vector<RDKit::Atom *> atomsNMV;
    GetAtomsWithNotMaxValence(mol, atomsNMV);

    std::vector<RDKit::Atom *> &atoms1 = atomsNMV;
    std::vector<RDKit::Atom *> &atoms2 = atomsNMV;

    for (int i = 0; i < atoms1.size(); ++i) {
        for (int j = 0; j < atoms2.size(); ++j) {
            if (atoms1[i]->getIdx() != atoms2[j]->getIdx()) {
                addBondCandidates.push_back(
                    std::make_pair(atoms1[i]->getIdx(), atoms2[j]->getIdx()));
            }
        }
    }

    std::vector<int> ixRemove;
    for (int i = 0; i < addBondCandidates.size(); ++i) {
        AtomIdx atomIdx1 = addBondCandidates[i].first;
        AtomIdx atomIdx2 = addBondCandidates[i].second;
        for (int j = i + 1; j < addBondCandidates.size(); ++j) {
            if ((addBondCandidates[j].first == atomIdx1 && addBondCandidates[j].second == atomIdx2) ||
                    (addBondCandidates[j].second == atomIdx1 && addBondCandidates[j].first == atomIdx2)) {
                ixRemove.push_back(j);
            }
        }
    }

    std::sort(ixRemove.begin(), ixRemove.end());
    for (int i = ixRemove.size() - 1; i >= 0; i--) {
        addBondCandidates.erase(addBondCandidates.begin() + ixRemove[i]);
    }
}

void MorphingData::InitMutateAtom()
{
    // Determine possible substitutions.
    RDKit::Atom *atom;
    RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
    RDKit::ROMol::AtomIterator iter;
    for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
        atom = *iter;
        std::vector<MolpherAtom> sub;
        for (MolpherAtomIdx idx = 0; idx < atoms.size(); ++idx) {
            if ((atom->getAtomicNum() != atoms[idx].atomicNum) &&
                    (atom->getExplicitValence() <= GetMaxBondsMod(atoms[idx]))) {
                sub.push_back(atoms[idx]);
            }
        }
        mutateAtomCandidates.push_back(sub);
    }
}

void MorphingData::InitRemoveAtom()
{
    // Find boundary atoms.
    RDKit::Atom *atom;
    RDKit::ROMol::AtomIterator iter;
    for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
        atom = *iter;
        if (RDKit::queryAtomHeavyAtomDegree(atom) == 1) {
            removeAtomCandidates.push_back(atom->getIdx());
        }
    }
}

void MorphingData::InitRemoveBond()
{
    // Find high-order or ring bonds.
    RDKit::Bond *bond;
    RDKit::ROMol::BondIterator iter;
    for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
        bond = *iter;
        if ((RDKit::queryBondOrder(bond) > 1) || RDKit::queryIsBondInRing(bond)) {
            removeBondCandidates.push_back(bond->getIdx());
        }
    }
}

void MorphingData::InitInterlayAtom()
{
    RDKit::Bond *bond;
    RDKit::ROMol::BondIterator iter;
    for (MolpherAtomIdx idx = 0; idx < atoms.size(); ++idx) {
        for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
            bond = *iter;
            // the aromatic bonds are reduced to single and double bonds
            if ((RDKit::queryBondOrder(bond) * 2) <= GetMaxBondsMod(atoms[idx])) {
                interlayAtomCandidates[idx].push_back(bond->getIdx());
            }
        }
    }
}

void MorphingData::InitBondReroute()
{
    RDKit::Bond *bond;
    RDKit::Atom *atom0, *atom1, *bondAtoms[2];
    std::vector<RDKit::Atom *> candidates[2];
    std::queue<RDKit::Atom *> q;

    RDKit::ROMol::BondIterator iter;
    // for each bond in the molecule
    for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
        // for begin and end atom of the selected bond
        bond = *iter;
        for (int i = 0; i < 2; ++i) {
            candidates[i].clear();
            // set for which one re-routing candidates will be identified
            if (i == 0) {
                bondAtoms[0] = bond->getBeginAtom();
                bondAtoms[1] = bond->getEndAtom();
            } else {
                bondAtoms[0] = bond->getEndAtom();
                bondAtoms[1] = bond->getBeginAtom();
            }

            // q is used for breadth-first-search
            q.push(bondAtoms[1]);

            while (!q.empty()) {
                atom0 = q.front();
                q.pop();

                RDKit::ROMol::ADJ_ITER beg, end;
                boost::tie(beg, end) = mol.getAtomNeighbors(atom0);
                while (beg != end) {
                    atom1 = mol[*beg++].get();
                    if ((atom0->getIdx() == bondAtoms[1]->getIdx()) &&
                            (atom1->getIdx() == bondAtoms[0]->getIdx())) {
                        // the same bond as we calculate with
                        continue;
                    }
                    if (atom1->getIdx() == bondAtoms[1]->getIdx()) {
                        // if we returned to the bond being rerouted
                        continue;
                    }
                    bool alreadyCandidate = false;
                    for (int j = 0; j < candidates[i].size(); ++j) {
                        if (candidates[i][j]->getIdx() == atom1->getIdx()) {
                            alreadyCandidate = true;
                            break;
                        }
                    }
                    if (!alreadyCandidate) {
                        if (atom1->getIdx() != bondAtoms[0]->getIdx()) {
                            candidates[i].push_back(atom1);
                        }
                        q.push(atom1);
                    }
                }
            }

            // remove candidates which do not have enough free valence or the
            // bond already exists
            for (int j = candidates[i].size() - 1; j >= 0 ; --j) {
                int maxBond = GetMaxBondsMod(*candidates[i][j]);
                int cntAlreadyBonded = candidates[i][j]->getExplicitValence();
                int bondOrder = RDKit::queryBondOrder(bond);

                // there is already a bond between atoms
                bool existsBond = mol.getBondBetweenAtoms(
                    bondAtoms[0]->getIdx(), candidates[i][j]->getIdx());

                if (((cntAlreadyBonded + bondOrder) > maxBond) || existsBond) {
                    candidates[i].erase(candidates[i].begin() + j);
                }
            }
        }
        if ((candidates[0].size() == 0) && (candidates[1].size() == 0)) {
            continue;
        }

        RerouteCandidates rc;
        rc.bondIdx = bond->getIdx();

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < candidates[i].size(); ++j) {
                rc.candidates[i].push_back(candidates[i][j]->getIdx());
            }
        }
        bondRerouteCandidates.push_back(rc);
    }
}

void MorphingData::InitBondContraction()
{
    RDKit::Bond *bond;
    RDKit::ROMol::BondIterator iter;
    for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
        bond = *iter;
        RDKit::Atom *atomToStay = bond->getBeginAtom();
        RDKit::Atom *atomToRemove = bond->getEndAtom();

        if ((atomToStay->getAtomicNum() == atomToRemove->getAtomicNum()) &&
                (RDKit::queryAtomHeavyAtomDegree(atomToStay) > 1) &&
                (RDKit::queryAtomHeavyAtomDegree(atomToRemove) > 1) &&
                ((atomToStay->getExplicitValence() +
                    atomToRemove->getExplicitValence() -
                    (2 * RDKit::queryBondOrder(bond))) <=
                        GetMaxBondsMod(*atomToStay))) {
            bondContractionCandidates.push_back(bond->getIdx());
        }
    }
}

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
#include <map>

#include <GraphMol/GraphMol.h>

#include "core/misc/global_types.h"
#include "selectors/chemoper_selectors.h"

class MorphingData
{
public:
    MorphingData(
        RDKit::ROMol &molecule,
        RDKit::ROMol &target,
        const std::vector<ChemOperSelector> &operators
        );
    MorphingData(
            RDKit::ROMol &molecule,
            const std::vector<ChemOperSelector> &operators,
            const std::set<int>& fixed_atom_ids
    );
    ~MorphingData();

protected:
    void InitAddAtom();
    void InitAddBond();
    void InitMutateAtom();
    void InitRemoveAtom();
    void InitRemoveBond();
    void InitInterlayAtom();
    void InitBondReroute();
    void InitBondContraction();

public:
    RDKit::ROMol &mol;
    std::vector<MolpherAtom> atoms;
    std::vector<ChemOperSelector> operators;
    std::set<int> fixed_atom_ids;
    tbb::concurrent_hash_map<std::string, int> removed_atoms;

    typedef struct {
        BondIdx bondIdx;
        std::vector<AtomIdx> candidates[2];
    } RerouteCandidates;

    std::vector<AtomIdx> addAtomCandidates;
    std::vector<AtomIdx> removeAtomCandidates;
    std::vector<std::pair<AtomIdx, AtomIdx> > addBondCandidates;
    std::vector<BondIdx> removeBondCandidates;
    std::vector<std::vector<MolpherAtom> > mutateAtomCandidates;
    std::map<AtomIdx, std::vector<BondIdx> > interlayAtomCandidates;
    std::vector<RerouteCandidates> bondRerouteCandidates;
    std::vector<BondIdx> bondContractionCandidates;
};

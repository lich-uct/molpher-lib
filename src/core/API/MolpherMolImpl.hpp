/*
 Copyright (c) 2016 Martin Šícho

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

#ifndef MOLPHERMOLIMPL_HPP
#define	MOLPHERMOLIMPL_HPP

#include <memory>
#include <vector>

#include "core/data_structs/MolpherMolData.hpp"
#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"

class MolpherMol::MolpherMolImpl {
    
    friend class MolpherMol; // TODO: this should not be necassary -> change architecture
    
private:
    std::shared_ptr<ExplorationTree> tree;
    MolpherMolData data;
    std::unique_ptr<RDKit::ROMol> rd_mol;
	std::set<int> fixed_atoms;

public:
    MolpherMolImpl(const std::string& string_repr);
    MolpherMolImpl(const MolpherMolData& data);
    MolpherMolImpl(const MolpherMolImpl& other);
	MolpherMolImpl(
			const RDKit::RWMol& rd_mol
			, const std::string& formula
			, const std::string& parentSmile
			, const unsigned& oper
			, const double& dist
			, const double& distToClosestDecoy
			, const double& weight
			, const double& sascore
			, const std::set<int>& fixed_atoms
	);
    MolpherMolImpl();
    
    std::unique_ptr<MolpherMolImpl> copy() const;

    void initialize(const std::string &string_repr);
    void initialize(std::unique_ptr<RDKit::RWMol> mol);

	std::unique_ptr<ConcurrentMolVector> morph(
			const std::vector<ChemOperSelector>& operators
			, int cntMorphs
			, int threadCnt
			, FingerprintSelector fingerprintSelector
			, SimCoeffSelector simCoeffSelector
			, const MolpherMol& target
	);
	std::unique_ptr<ConcurrentMolVector> morph(
			const std::vector<ChemOperSelector>& operators
			, int cntMorphs
			, int threadCnt
	);
};

#endif	/* MOLPHERMOLIMPL_HPP */


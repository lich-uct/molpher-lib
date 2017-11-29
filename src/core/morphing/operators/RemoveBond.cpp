//
// Created by sichom on 11/26/17.
//

#include "RemoveBondImpl.hpp"

#include "core/chem/ChemicalAuxiliary.h"
#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

RemoveBond::RemoveBond() :
		MorphingOperator()
		, pimpl(new RemoveBondImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

RemoveBond::RemoveBondImpl::RemoveBondImpl() :
		MorphingOperatorImpl()
		, original_rdkit(nullptr)
{
	// no action
}

void RemoveBond::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> RemoveBond::morph() {
	return pimpl->morph();
}

const std::vector<std::pair<unsigned int, unsigned int>> &RemoveBond::getOpenBonds() {
	return pimpl->getOpenBonds();
}

void RemoveBond::RemoveBondImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		open_bonds.clear();
		open_bonds_rd.clear();

		// Find high-order or ring bonds.
		RDKit::Bond *bond;
		RDKit::ROMol::BondIterator iter;
		for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
			bond = *iter;
			AtomIdx begin_atm = bond->getBeginAtom()->getIdx();
			AtomIdx end_atm = bond->getEndAtom()->getIdx();

			if (
					(original->getAtom(begin_atm)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS_AND_BONDS)
				|| (original->getAtom(end_atm)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS_AND_BONDS)
					) {
				continue;
			}

			if ((RDKit::queryBondOrder(bond) > 1) || RDKit::queryIsBondInRing(bond)) {
				open_bonds_rd.push_back(bond->getIdx());
				open_bonds.push_back(
						std::make_pair(
								begin_atm
								, end_atm
						)
				);
			}
		}

		assert(open_bonds_rd.size() == open_bonds.size());
	} else {
		std::runtime_error("Invalid reference to original molecule.");
	}
}

std::shared_ptr<MolpherMol> RemoveBond::RemoveBondImpl::morph() {
	RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

	if (open_bonds_rd.size() == 0) {
		delete newMol;
		SynchCerr("No open atom pairs for bond addition. Skipping...");
		return nullptr;
	}

	int randPos = SynchRand::GetRandomNumber(open_bonds_rd.size() - 1);
	RDKit::Bond *bond = newMol->getBondWithIdx(open_bonds_rd[randPos]);

	if (RDKit::queryIsBondInRing(bond)) {
		newMol->removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
	} else {
		DecreaseBondOrder(*bond);
	}

	std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
	writeOriginalLockInfo(ret);

	return ret;
}

const std::vector<std::pair<AtomIdx, AtomIdx>>& RemoveBond::RemoveBondImpl::getOpenBonds() {
	return open_bonds;
}

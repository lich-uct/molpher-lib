//
// Created by sichom on 11/22/17.
//

#include "AddBondImpl.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <core/chem/ChemicalAuxiliary.h>
#include <core/misc/inout.h>
#include <core/misc/SynchRand.h>

AddBond::AddBond() :
		MorphingOperator()
		, pimpl(new AddBondImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

AddBond::AddBondImpl::AddBondImpl() :
		MorphingOperatorImpl()
		, original_rdkit(nullptr)
{
	// no action
}

void AddBond::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> AddBond::morph() {
	return pimpl->morph();
}

const std::vector<std::pair<unsigned int, unsigned int>> &AddBond::getOpenBonds() {
	return pimpl->getOpenBonds();
}

void AddBond::AddBondImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		open_bonds.clear();

		std::vector<RDKit::Atom *> atomsNMV;
		RDKit::Atom *atom;
		RDKit::ROMol::AtomIterator iter;
		for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
			atom = *iter;
			if (original->getAtom(atom->getIdx())->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS_AND_BONDS) {
				// pass if the bonds should not be modified
				continue;
			}

			if (CntFreeBonds(*atom)) {
				atomsNMV.push_back(atom);
			}
		}

		std::vector<RDKit::Atom *> &atoms1 = atomsNMV;
		std::vector<RDKit::Atom *> &atoms2 = atomsNMV;
		for (int i = 0; i < atoms1.size(); ++i) {
			for (int j = 0; j < atoms2.size(); ++j) {
				if (atoms1[i]->getIdx() != atoms2[j]->getIdx()) {
					open_bonds.push_back(
							std::make_pair(atoms1[i]->getIdx(), atoms2[j]->getIdx()));
				}
			}
		}

		std::vector<int> ixRemove;
		for (int i = 0; i < open_bonds.size(); ++i) {
			AtomIdx atomIdx1 = open_bonds[i].first;
			AtomIdx atomIdx2 = open_bonds[i].second;
			for (int j = i + 1; j < open_bonds.size(); ++j) {
				if ((open_bonds[j].first == atomIdx1 && open_bonds[j].second == atomIdx2) ||
					(open_bonds[j].second == atomIdx1 && open_bonds[j].first == atomIdx2)) {
					ixRemove.push_back(j);
				}
			}
		}

		std::sort(ixRemove.begin(), ixRemove.end());
		for (int i = ixRemove.size() - 1; i >= 0; i--) {
			open_bonds.erase(open_bonds.begin() + ixRemove[i]);
		}
	} else {
		std::runtime_error("Invalid reference for original molecule.");
	}
}

std::shared_ptr<MolpherMol> AddBond::AddBondImpl::morph() {
	RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

	if (open_bonds.size() == 0) {
		delete newMol;
		SynchCerr("No open atom pairs for bond addition. Skipping...");
		return nullptr;
	}

	int randPos = SynchRand::GetRandomNumber(open_bonds.size() - 1);
	AtomIdx idx1 = open_bonds[randPos].first;
	AtomIdx idx2 = open_bonds[randPos].second;
	RDKit::Bond *bond = newMol->getBondBetweenAtoms(idx1, idx2);
	if (!bond) {
		newMol->addBond(idx1, idx2, RDKit::Bond::SINGLE);
	} else {
		IncreaseBondOrder(*bond);
	}

	std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
	writeOriginalLockInfo(ret);

	return ret;
}

const std::vector<std::pair<AtomIdx, AtomIdx>>& AddBond::AddBondImpl::getOpenBonds() {
	return open_bonds;
}

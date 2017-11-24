//
// Created by sichom on 7/25/17.
//

#include "AddAtomImpl.hpp"
#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

AddAtom::AddAtom() :
MorphingOperator()
, pimpl(new AddAtomImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

void AddAtom::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> AddAtom::morph() {
	return pimpl->morph();
}

const std::vector<unsigned int>& AddAtom::getOpenIndices() {
	return pimpl->getOpenIndices();
}

std::vector<std::shared_ptr<MolpherAtom>> AddAtom::getOpenAtoms() {
	return pimpl->getOpenAtoms();
}

AddAtom::AddAtom(const AtomLibrary &atom_library) :
MorphingOperator()
, pimpl(new AddAtomImpl(atom_library))
{
	setMorphingOperatorPimpl(pimpl);
}

AddAtom::AddAtomImpl::AddAtomImpl() :
MorphingOperatorImpl()
, original_rdkit(nullptr)
, atom_library(AtomLibrary::getDefaultLibrary())
{
	// no action
}

void AddAtom::AddAtomImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	// FIXME: check for validity of the pointer
	original = mol_orig;
	original_rdkit.reset(original->asRDMol());
	RDKit::ROMol& mol = *original_rdkit;
	open_atoms.clear();

	RDKit::Atom *atom;
	RDKit::ROMol::AtomIterator iter;
	for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
		atom = *iter;
		int free_bonds = atom->getImplicitValence();

		if (free_bonds >= 1 && !(MolpherAtom::NO_ADDITION & original->getAtom(atom->getIdx())->getLockingMask())) {
			open_atoms.push_back(atom->getIdx());
		}
	}
}

std::shared_ptr<MolpherMol> AddAtom::AddAtomImpl::morph() {
	if (original_rdkit) {
		RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

		RDKit::Atom* atom = atom_library.getRandomAtom().asRDAtom();

		if (open_atoms.size() == 0) {
			delete newMol;
			SynchCerr("No open atoms for addition. Skipping...");
			return nullptr;
		}

		AtomIdx bindingAtomIdx = open_atoms[SynchRand::GetRandomNumber(open_atoms.size() - 1)];
		RDKit::Atom *bindingAtom = newMol->getAtomWithIdx(bindingAtomIdx);

		AtomIdx newAtomIdx = newMol->addAtom(atom);

		if (HasNonSingleBond(*bindingAtom)
			&& !(original->getAtom(bindingAtomIdx)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS_AND_BONDS)
			&& (SynchRand::GetRandomNumber(0, 1) > 0)
				) {
			RDKit::Bond *bond = GetRandomNonSingleBond(*bindingAtom);
			DecreaseBondOrder(*bond);
		}

		newMol->addBond(bindingAtom->getIdx(), newAtomIdx, RDKit::Bond::SINGLE);
		bindingAtom->calcExplicitValence();
		bindingAtom->calcImplicitValence();

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret);

		delete atom;

		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original molecule to morph first.");
	}
}

const std::vector<unsigned int> &AddAtom::AddAtomImpl::getOpenIndices() {
	return open_atoms;
}

std::vector<std::shared_ptr<MolpherAtom>> AddAtom::AddAtomImpl::getOpenAtoms() {
	std::vector<std::shared_ptr<MolpherAtom>> ret;
	for (auto idx : open_atoms) {
		ret.push_back(original->getAtom(idx));
	}
	return ret;
}

AddAtom::AddAtomImpl::AddAtomImpl(const AtomLibrary &atom_library) :
AddAtomImpl()
{
	this->atom_library = atom_library;
}

//
// Created by sichom on 11/28/17.
//

#include "InterlayAtomImpl.hpp"
#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

InterlayAtom::InterlayAtom() :
		MorphingOperator()
		, pimpl(new InterlayAtomImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

void InterlayAtom::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> InterlayAtom::morph() {
	return pimpl->morph();
}

InterlayAtom::InterlayAtom(const AtomLibrary &atom_library) :
		MorphingOperator()
		, pimpl(new InterlayAtomImpl(atom_library))
{
	setMorphingOperatorPimpl(pimpl);
}

InterlayAtom::InterlayAtomImpl::InterlayAtomImpl() :
		MorphingOperatorImpl()
		, atom_library(AtomLibrary::getDefaultLibrary())
{
	// no action
}

void InterlayAtom::InterlayAtomImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		MorphingOperatorImpl::setOriginal(mol_orig);

		RDKit::ROMol& mol = *original_rdkit;
		interlay_candidates.clear();
		atoms = atom_library.getAtoms();

		RDKit::Bond *bond;
		RDKit::ROMol::BondIterator iter;
		for (AtomIdx idx = 0; idx < atoms.size(); ++idx) {
			for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
				bond = *iter;

				AtomIdx begin_idx = bond->getBeginAtomIdx();
				AtomIdx end_idx = bond->getEndAtomIdx();
				if (
					(original->getAtom(begin_idx)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS)
					|| (original->getAtom(begin_idx)->getLockingMask() & MolpherAtom::KEEP_BONDS)
					|| (original->getAtom(end_idx)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS)
					|| (original->getAtom(end_idx)->getLockingMask() & MolpherAtom::KEEP_BONDS)
				) {
					continue;
				}

				// the aromatic bonds are reduced to single and double bonds
				if ((RDKit::queryBondOrder(bond) * 2) <= GetMaxBondsMod(*atoms[idx])) {
					interlay_candidates[idx].push_back(bond->getIdx());
				}
			}
		}
	} else {
		throw std::runtime_error("Invalid reference to original molecule.");
	}
}

std::shared_ptr<MolpherMol> InterlayAtom::InterlayAtomImpl::morph() {
	if (original_rdkit) {
		auto *newMol(new RDKit::RWMol(*original_rdkit));

		RDKit::Atom atom;
		AtomIdx idx = GetRandomAtom(atoms, atom);

		if (interlay_candidates.find(idx) == interlay_candidates.end()) {
			delete newMol;
//			SynchCerr("Given bond cannot be interlayed with the selected atom (" + atom.getSymbol() + "). Skipping: " + original->getSMILES());
			return nullptr;
		}
		if (interlay_candidates[idx].empty()) {
			delete newMol;
//			SynchCerr("No bond to interlay.  Skipping: " + original->getSMILES());
			return nullptr;
		}


		int randPos = SynchRand::GetRandomNumber(interlay_candidates[idx].size() - 1);
		RDKit::Bond *bond = newMol->getBondWithIdx(interlay_candidates[idx][randPos]);

		AtomIdx beg = bond->getBeginAtomIdx();
		AtomIdx end = bond->getEndAtomIdx();
		RDKit::Bond::BondType bt = bond->getBondType(); // FIXME: this could be a problem in rings (two double bonds next to each other for example)

		AtomIdx newAtomIdx = newMol->addAtom(&atom); // atom is copied

		newMol->removeBond(beg, end);
		newMol->addBond(beg, newAtomIdx, bt);
		newMol->addBond(newAtomIdx, end, bt);

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret);
		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original molecule to morph first.");
	}
}

InterlayAtom::InterlayAtomImpl::InterlayAtomImpl(const AtomLibrary &atom_library) :
		InterlayAtomImpl()
{
	this->atom_library = atom_library;
}

std::string InterlayAtom::getName() const {
	return ChemOperLongDesc(OP_INTERLAY_ATOM);
}
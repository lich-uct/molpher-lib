//
// Created by sichom on 8/28/17.
//

#include <GraphMol/ROMol.h>
#include <GraphMol/QueryOps.h>
#include <core/misc/inout.h>
#include <core/misc/SynchRand.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "RemoveAtomImpl.hpp"

RemoveAtom::RemoveAtom() :
MorphingOperator()
, pimpl(new RemoveAtomImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

void RemoveAtom::setOriginal(std::shared_ptr<MolpherMol> mol)
{
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> RemoveAtom::morph() {
	return pimpl->morph();
}

const std::vector<unsigned int> &RemoveAtom::getMarkedIndices() {
	return pimpl->getMarkedIndices();
}

std::vector<std::shared_ptr<MolpherAtom>> RemoveAtom::getMarkedAtoms() {
	return pimpl->getMarkedAtoms();
}

// implementation

RemoveAtom::RemoveAtomImpl::RemoveAtomImpl() :
MorphingOperatorImpl()
, original_rdkit(nullptr)
{
	// no action
}

void RemoveAtom::RemoveAtomImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		marked_atoms.clear();

		// Find boundary atoms.
		RDKit::Atom *atom;
		RDKit::Atom *neighbor;
		RDKit::ROMol::AtomIterator iter;
		for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
			atom = *iter;

//		bool kept_by_neighbor = false;
//		RDKit::ROMol::ADJ_ITER beg, end;
//		boost::tie(beg, end) = mol.getAtomNeighbors(atom);
//		while (beg != end) {
//			neighbor = mol[*beg].get();
//			int locking_mask = original->getAtom(neighbor->getIdx())->getLockingMask();
//			if (MolpherAtom::KEEP_NEIGHBORS & locking_mask) {
//				kept_by_neighbor = true;
//				break;
//			}
//			++beg;
//		}

			if (RDKit::queryAtomHeavyAtomDegree(atom) == 1
				//			&& !kept_by_neighbor
				&& !((MolpherAtom::NO_REMOVAL) & original->getAtom(atom->getIdx())->getLockingMask())) {
				marked_atoms.push_back(atom->getIdx());
			}
		}
	} else {
		std::runtime_error("Invalid reference for original molecule.");
	}
}

std::shared_ptr<MolpherMol> RemoveAtom::RemoveAtomImpl::morph() {
	if (original_rdkit) {
		RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

		if (marked_atoms.size() == 0) {
			delete newMol;
			SynchCerr("No atoms marked for removal. Skipping...");
			return nullptr;
		}

		// remove a random boundary atom
		int randPos = SynchRand::GetRandomNumber(marked_atoms.size() - 1);
		AtomIdx removedIdx = marked_atoms[randPos];
		newMol->removeAtom(removedIdx);

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret, removedIdx);

		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original molecule to morph first.");
	}
}

const std::vector<unsigned int> &RemoveAtom::RemoveAtomImpl::getMarkedIndices() {
	return marked_atoms;
}

std::vector<std::shared_ptr<MolpherAtom>> RemoveAtom::RemoveAtomImpl::getMarkedAtoms() {
	std::vector<std::shared_ptr<MolpherAtom>> ret;
	for (auto idx : marked_atoms) {
		ret.push_back(original->getAtom(idx));
	}
	return ret;
}

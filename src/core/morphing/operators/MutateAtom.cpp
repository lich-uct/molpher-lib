//
// Created by sichom on 11/27/17.
//

#include "MutateAtomImpl.hpp"
#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>

MutateAtom::MutateAtom() :
		MorphingOperator()
		, pimpl(new MutateAtomImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

void MutateAtom::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> MutateAtom::morph() {
	return pimpl->morph();
}

MutateAtom::MutateAtom(const AtomLibrary &atom_library) :
		MorphingOperator()
		, pimpl(new MutateAtomImpl(atom_library))
{
	setMorphingOperatorPimpl(pimpl);
}

std::string MutateAtom::getName() const {
	return "Mutate Atom";
}

MutateAtom::MutateAtomImpl::MutateAtomImpl() :
		MorphingOperatorImpl()
		, original_rdkit(nullptr)
		, atom_library(AtomLibrary::getDefaultLibrary())
{
	// no action
}

void MutateAtom::MutateAtomImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		replacements.clear();

		RDKit::Atom *atom_rd;
		RDKit::ROMol::AtomIterator iter;
		auto atoms = atom_library.getAtoms();
		for (iter = mol.beginAtoms(); iter != mol.endAtoms(); iter++) {
			atom_rd = *iter;
			MolpherAtom& atom = *original->getAtom(atom_rd->getIdx());
			std::vector<MolpherAtom> atom_replacements;
//			SynchCout(atom.getSymbol());
//			std::cout << atom_rd->getExplicitValence() << std::endl;

			if (atom.getLockingMask() & MolpherAtom::NO_MUTATION) {
				replacements.push_back(atom_replacements);
				continue;
			}

//			RDKit::ROMol::OEDGE_ITER beg,end;
//			boost::tie(beg,end) = mol.getAtomBonds(atom_rd);
//			bool locked_by_neighbor = false;
//			while(beg!=end){
//				const RDKit::BOND_SPTR bond = (mol)[*beg];
//
//				AtomIdx other_idx = bond->getOtherAtomIdx(atom_rd->getIdx());
//				if (
//						(original->getAtom(other_idx)->getLockingMask() & MolpherAtom::KEEP_NEIGHBORS)
//						|| (original->getAtom(other_idx)->getLockingMask() & MolpherAtom::KEEP_BONDS)
//						) {
//					locked_by_neighbor = true;
//					break;
//				}
//
//				++beg;
//			}
//			if (locked_by_neighbor) {
//				replacements.push_back(atom_replacements);
//				continue;
//			}

			for (AtomIdx idx = 0; idx < atoms.size(); ++idx) {
				if ((atom_rd->getAtomicNum() != atoms[idx]->getAtomicNum()) &&
					(atom_rd->getExplicitValence() <= GetMaxBondsMod(*atoms[idx]))) {
					// TODO: is this what we want? -- stereochemistry plays a role here -> if stereochemistry is defined, an implicit hydrogen is considered explicit
					atom_replacements.push_back(*atoms[idx]);
				}
			}

			replacements.push_back(atom_replacements);
		}

		assert(replacements.size() == original->getAtomCount());
	} else {
		throw std::runtime_error("Invalid reference to original molecule.");
	}
}

std::shared_ptr<MolpherMol> MutateAtom::MutateAtomImpl::morph() {
	if (original_rdkit) {
		RDKit::RWMol *newMol = original->asRDMol();

		int randPos = SynchRand::GetRandomNumber(newMol->getNumAtoms() - 1);

		if(replacements[randPos].empty()) {
			delete newMol;
//			SynchCerr("Given atom cannot be mutated.  Skipping: " + original->getSMILES());
			return nullptr;
		}

		RDKit::Atom* atom = GetRandomAtom(atom_library);
		bool validReplacement = false;
		for (auto& atm : replacements[randPos]) {
			if ((atm.getAtomicNum() == atom->getAtomicNum()) && (atm.getFormalCharge() == atom->getFormalCharge())) {
				validReplacement = true;
//				SynchCout("Found valid replacement: " + atm.getSymbol());
				break;
			}
		}
		if (!validReplacement) {
//			SynchCerr("Generated replacement (" + atom->getSymbol() + ") cannot be applied.  Skipping: " + original->getSMILES());
			delete newMol;
			delete atom;
			return nullptr;
		}

		newMol->replaceAtom(randPos, atom); // atom is copied
		delete atom;

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret);
		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original molecule to morph first.");
	}
}

MutateAtom::MutateAtomImpl::MutateAtomImpl(const AtomLibrary &atom_library) :
		MutateAtomImpl()
{
	this->atom_library = atom_library;
}
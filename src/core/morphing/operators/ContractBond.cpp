//
// Created by sichom on 11/29/17.
//

#include "ContractBondImpl.hpp"

#include "core/chem/ChemicalAuxiliary.h"
#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"

ContractBond::ContractBond() :
		MorphingOperator()
		, pimpl(new ContractBondImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

ContractBond::ContractBondImpl::ContractBondImpl() :
		MorphingOperatorImpl()
		, original_rdkit(nullptr)
{
	// no action
}

void ContractBond::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> ContractBond::morph() {
	return pimpl->morph();
}

const std::vector<std::pair<unsigned int, unsigned int>> &ContractBond::getOpenBonds() {
	return pimpl->getOpenBonds();
}

std::string ContractBond::getName() const {
	return ChemOperLongDesc(OP_BOND_CONTRACTION);
}

void ContractBond::ContractBondImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		open_bonds.clear();
		open_bonds_rd.clear();

		RDKit::ROMol::BondIterator iter;
		for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
			RDKit::Bond *bond = *iter;

			tryToOpenBond(bond);
			tryToOpenBond(bond, true);
		}

		assert(open_bonds_rd.size() == open_bonds.size());
	} else {
		throw std::runtime_error("Invalid reference to original molecule.");
	}
}

void ContractBond::ContractBondImpl::tryToOpenBond(RDKit::Bond *bond, bool swap) {
	RDKit::Atom *atomToStay = nullptr;
	RDKit::Atom *atomToRemove = nullptr;
	if (!swap) {
		atomToStay = bond->getBeginAtom();
		atomToRemove = bond->getEndAtom();
	} else {
		atomToStay = bond->getEndAtom();
		atomToRemove = bond->getBeginAtom();
	}
	AtomIdx atomToStayIdx = atomToStay->getIdx();
	AtomIdx atomToRemoveIdx = atomToRemove->getIdx();

	RDKit::ROMol& mol = *original_rdkit;
	if ((atomToStay->getAtomicNum() == atomToRemove->getAtomicNum()) &&
		(RDKit::queryAtomHeavyAtomDegree(atomToStay) > 1) &&
		(RDKit::queryAtomHeavyAtomDegree(atomToRemove) > 1) &&
		((atomToStay->getExplicitValence() +
		  atomToRemove->getExplicitValence() -
		  (2 * RDKit::queryBondOrder(bond))) <=
		 GetMaxBondsMod(*atomToStay))) {

		bool no_removal = false;
		if (original->getAtom(atomToRemoveIdx)->getLockingMask() & MolpherAtom::NO_REMOVAL) {
			no_removal = true;
		}

//		bool kept_by_neighbor = false;
//		if (!no_removal) {
//			RDKit::Atom *neighbor;
//			RDKit::ROMol::ADJ_ITER beg, end;
//			boost::tie(beg, end) = mol.getAtomNeighbors(atomToRemove);
//			while (beg != end) {
//				neighbor = mol[*beg].get();
//				MolpherAtom::LockingMask locking_mask = (MolpherAtom::LockingMask) original->getAtom(neighbor->getIdx())->getLockingMask();
//				if (MolpherAtom::KEEP_NEIGHBORS & locking_mask) {
//					kept_by_neighbor = true;
//					break;
//				}
//				++beg;
//			}
//		}


		if (!no_removal) {
			open_bonds_rd.push_back(bond->getIdx());
			open_bonds.push_back(std::make_pair(atomToStayIdx, atomToRemoveIdx));
		}
	}
}

std::shared_ptr<MolpherMol> ContractBond::ContractBondImpl::morph() {
	if (original_rdkit) {
		RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

		if (open_bonds_rd.size() == 0) {
			delete newMol;
//			SynchCerr("No bonds open for contraction. Skipping: " + original->getSMILES());
			return nullptr;
		}

		int randPos = SynchRand::GetRandomNumber(open_bonds.size() - 1);
		RDKit::Bond *bond = newMol->getBondWithIdx(open_bonds_rd[randPos]);
		RDKit::Atom *atomToRemove = nullptr;
		RDKit::Atom *atomToStay = nullptr;
		AtomIdx atomToStayIdx = open_bonds[randPos].first;
		AtomIdx atomToRemoveIdx = open_bonds[randPos].second;
		if (bond->getBeginAtomIdx() == atomToStayIdx) {
			atomToStay = bond->getBeginAtom();
			atomToRemove = bond->getEndAtom();
		} else {
			atomToStay = bond->getEndAtom();
			atomToRemove = bond->getBeginAtom();
		}


		std::vector<RDKit::Bond *> bondsToRemove;
		RDKit::Bond *bondToChange;
		RDKit::ROMol::OEDGE_ITER beg, end;
		boost::tie(beg, end) = newMol->getAtomBonds(atomToRemove);
		while (beg != end) {
			bondToChange = (*original_rdkit)[*beg++].get();
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


		for (int i = 0; i < bondsToRemove.size(); ++i) {
			newMol->removeBond(
					bondsToRemove[i]->getBeginAtomIdx()
					, bondsToRemove[i]->getEndAtomIdx()
			);
		}
		newMol->removeAtom(atomToRemove);

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret, atomToRemoveIdx);

		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original structure first.");
	}
}

const std::vector<std::pair<AtomIdx, AtomIdx>>& ContractBond::ContractBondImpl::getOpenBonds() {
	return open_bonds;
}
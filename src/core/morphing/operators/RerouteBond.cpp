//
// Created by sichom on 1/8/18.
//

//
// Created by sichom on 11/22/17.
//

#include "RerouteBondImpl.hpp"

#include <queue>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <core/chem/ChemicalAuxiliary.h>
#include <core/misc/inout.h>
#include <core/misc/SynchRand.h>

RerouteBond::RerouteBond() :
		MorphingOperator()
		, pimpl(new RerouteBondImpl())
{
	setMorphingOperatorPimpl(pimpl);
}

RerouteBond::RerouteBondImpl::RerouteBondImpl() :
		MorphingOperatorImpl()
		, original_rdkit(nullptr)
{
	// no action
}

void RerouteBond::setOriginal(std::shared_ptr<MolpherMol> mol) {
	pimpl->setOriginal(mol);
}

std::shared_ptr<MolpherMol> RerouteBond::morph() {
	return pimpl->morph();
}

std::string RerouteBond::getName() const {
	return ChemOperLongDesc(OP_BOND_REROUTE);
}

void RerouteBond::RerouteBondImpl::setOriginal(std::shared_ptr<MolpherMol> mol_orig) {
	if (mol_orig) {
		original = mol_orig;
		original_rdkit.reset(original->asRDMol());
		RDKit::ROMol& mol = *original_rdkit;
		candidates.clear();

		RDKit::Bond *bond;
		RDKit::Atom *atom0, *atom1, *bondAtoms[2];
		std::vector<RDKit::Atom *> candidates[2];
		std::queue<RDKit::Atom *> q;

		RDKit::ROMol::BondIterator iter;
		// for each bond in the molecule
		for (iter = mol.beginBonds(); iter != mol.endBonds(); iter++) {
			// for begin and end atom of the selected bond
			bond = *iter;
			int begin_lock = original->getAtom(bond->getBeginAtom()->getIdx())->getLockingMask();
			int end_lock = original->getAtom(bond->getEndAtom()->getIdx())->getLockingMask();
			int mask = MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::KEEP_BONDS;
			if ((begin_lock & mask) || (end_lock & mask)) {
				continue;
			}

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

						// do not add candidates with locked atom additions (the bond should not be rerouted to those because that would create an additional bond on the new neighbor (new bond end atom), which is not allowed)
						if (original->getAtom(atom1->getIdx())->getLockingMask() & MolpherAtom::NO_ADDITION) {
							continue;
						}

						// push to candidates if not there already
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

			Candidates rc;
			rc.bondIdx = bond->getIdx();

			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < candidates[i].size(); ++j) {
					rc.candidates[i].push_back(candidates[i][j]->getIdx());
				}
			}
			this->candidates.push_back(rc);
		}
	} else {
		throw std::runtime_error("Invalid reference to original molecule.");
	}
}

std::shared_ptr<MolpherMol> RerouteBond::RerouteBondImpl::morph() {
	if (original_rdkit) {
		RDKit::RWMol *newMol = new RDKit::RWMol(*original_rdkit);

		if (candidates.size() == 0) {
			delete newMol;
//			SynchCerr("No candidates for bond reroute identified. Skipping: " + original->getSMILES());
			return nullptr;
		}
		int randPos = SynchRand::GetRandomNumber(candidates.size() - 1);
		RDKit::Bond *bond = newMol->getBondWithIdx(candidates[randPos].bondIdx);
		int randSide = SynchRand::GetRandomNumber(0, 1);
		if (candidates[randPos].candidates[randSide].size() == 0) {
			randSide = 1 - randSide;
		}

		int randPos2 =
				SynchRand::GetRandomNumber(
						candidates[randPos].candidates[randSide].size() - 1);

		AtomIdx beginAtomIdx =
				(randSide == 0) ? bond->getBeginAtomIdx() : bond->getEndAtomIdx();
		AtomIdx endAtomIdx =
				candidates[randPos].candidates[randSide][randPos2];
		AtomIdx remBeginAtom = bond->getBeginAtomIdx();
		AtomIdx remEndAtom = bond->getEndAtomIdx();

		newMol->addBond(beginAtomIdx, endAtomIdx, bond->getBondType());
		newMol->removeBond(remBeginAtom, remEndAtom);

		std::shared_ptr<MolpherMol> ret(new MolpherMol(newMol));
		writeOriginalLockInfo(ret);

		return ret;
	} else {
		throw std::runtime_error("No starting molecule set. Set the original structure first.");
	}
}

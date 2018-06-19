//
// Created by sichom on 7/25/17.
//

#include "MolpherAtomImpl.hpp"

#include "GraphMol/Atom.h"

const std::vector<MolpherAtom::LockingMask> MolpherAtom::atom_locks = {
UNLOCKED
, NO_MUTATION
, NO_ADDITION
, NO_REMOVAL
, KEEP_NEIGHBORS
, KEEP_NEIGHBORS_AND_BONDS
, KEEP_BONDS
, FULL_LOCK
};

const std::map<MolpherAtom::LockingMask, std::string> MolpherAtom::MolpherAtomImpl::locks_map = {
{MolpherAtom::UNLOCKED,"UNLOCKED"}
, {MolpherAtom::NO_ADDITION,"NO_ADDITION"}
, {MolpherAtom::NO_REMOVAL,"NO_REMOVAL"}
, {MolpherAtom::NO_MUTATION,"NO_MUTATION"}
, {MolpherAtom::KEEP_NEIGHBORS,"KEEP_NEIGHBORS"}
, {MolpherAtom::KEEP_NEIGHBORS_AND_BONDS,"KEEP_NEIGHBORS_AND_BONDS"}
, {MolpherAtom::KEEP_BONDS,"KEEP_BONDS"}
, {MolpherAtom::FULL_LOCK,"FULL_LOCK"}
};

MolpherAtom::MolpherAtomImpl::MolpherAtomImpl(RDKit::Atom *rd_atom)
:
atom_rd(rd_atom->copy())
, locking_mask(0)
{
	// no action
}

MolpherAtom::MolpherAtomImpl::MolpherAtomImpl(const std::string& symbol)
:
atom_rd(new RDKit::Atom(symbol))
, locking_mask(0)
{
	// no action
}

MolpherAtom::MolpherAtom(RDKit::Atom *rd_atom) :
pimpl(new MolpherAtomImpl(rd_atom))
{
	// no action
}

bool MolpherAtom::isLocked() const {
	return pimpl->locking_mask != 0;
}

int MolpherAtom::getLockingMask() const {
	return pimpl->locking_mask;
}

void MolpherAtom::setLockingMask(int mask) {
	if (mask & KEEP_NEIGHBORS_AND_BONDS) {
		mask = mask | KEEP_BONDS;
		mask = mask | KEEP_NEIGHBORS;
	}
	if (mask & (KEEP_BONDS | KEEP_NEIGHBORS)) {
		mask = mask | NO_REMOVAL; // TODO: does this make sense, should keeping neighbors/bonds imply not removing the atom itself?
	}
	pimpl->locking_mask = mask;
}

unsigned int MolpherAtom::getAtomicNum() const {
	return (unsigned int) pimpl->atom_rd->getAtomicNum();
}

double MolpherAtom::getMass() const {
	return pimpl->atom_rd->getMass();
}

int MolpherAtom::getFormalCharge() const {
	return pimpl->atom_rd->getFormalCharge();
}

MolpherAtom::MolpherAtom(const std::string& symbol) :
pimpl(new MolpherAtomImpl(symbol))
{
	// no action
}

MolpherAtom::MolpherAtom(const std::string& symbol, int charge) :
MolpherAtom(symbol)
{
	setFormalCharge(charge);
}

RDKit::Atom *MolpherAtom::asRDAtom() const {
	RDKit::Atom* new_atom = new RDKit::Atom(this->getAtomicNum());
	new_atom->setFormalCharge(this->getFormalCharge());
	return new_atom;
}

std::string MolpherAtom::getSymbol() const {
	return pimpl->atom_rd->getSymbol();
}

void MolpherAtom::setFormalCharge(int charge) {
	pimpl->atom_rd->setFormalCharge(charge);
}

MolpherAtom::MolpherAtom(const MolpherAtom &other)
: pimpl(new MolpherAtomImpl(other.getSymbol()))
{
	pimpl->atom_rd.reset(other.asRDAtom()); // FIXME: this is a result of a crappy hack -> atom_rd should only be initialized once
	this->setFormalCharge(other.getFormalCharge());
	this->setLockingMask(other.getLockingMask());
}

MolpherAtom::~MolpherAtom() = default;

MolpherAtom &MolpherAtom::operator=(const MolpherAtom& other) {
	pimpl = std::move(std::unique_ptr<MolpherAtomImpl>(new MolpherAtomImpl(other.asRDAtom())));
	this->setFormalCharge(other.getFormalCharge());
	this->setLockingMask(other.getLockingMask());
	return *this;
}

std::string MolpherAtom::lockToString(int lock) {
	// TODO: throw an exception if the lock is not found
	return MolpherAtomImpl::locks_map.find((LockingMask) lock)->second;
}

std::vector<std::string> MolpherAtom::lockingMaskToString(int mask) {
	std::vector<std::string> ret;
	for (LockingMask lock : atom_locks) {
		if (lock & mask && lock != FULL_LOCK) {
			ret.push_back(lockToString(lock));
		}
	}
	return ret;
}

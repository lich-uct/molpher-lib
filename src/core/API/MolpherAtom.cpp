//
// Created by sichom on 7/25/17.
//

#include "MolpherAtomImpl.hpp"

#include "GraphMol/Atom.h"

MolpherAtom::MolpherAtomImpl::MolpherAtomImpl(RDKit::Atom *rd_atom)
:
atomic_num((unsigned int) rd_atom->getAtomicNum())
, formal_charge(rd_atom->getFormalCharge())
, mass(rd_atom->getMass())
, symbol(rd_atom->getSymbol())
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
	pimpl->locking_mask = mask;
}

unsigned int MolpherAtom::getAtomicNum() const {
	return pimpl->atomic_num;
}

double MolpherAtom::getMass() const {
	return pimpl->mass;
}

int MolpherAtom::getFormalCharge() const {
	return pimpl->formal_charge;
}

MolpherAtom::MolpherAtom(std::string symbol) :
MolpherAtom(new RDKit::Atom(symbol))
{
	// no action
}

MolpherAtom::MolpherAtom(std::string symbol, int charge) :
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
	return pimpl->symbol;
}

void MolpherAtom::setFormalCharge(int charge) {
	pimpl->formal_charge = charge;
}

MolpherAtom::MolpherAtom(const MolpherAtom &other)
: pimpl(new MolpherAtomImpl(other.asRDAtom()))
{
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

//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_MOLPHERATOM_HPP
#define MOLPHER_LIB_MOLPHERATOM_HPP

#include <memory>
#include <vector>

namespace RDKit {
	class Atom;
}

class MolpherAtom {
private:
	class MolpherAtomImpl;
	std::unique_ptr<MolpherAtomImpl> pimpl;

public:
	enum LockingMask {
		UNLOCKED = 0 // everything goes (the default option)
		, NO_MUTATION = 1<<1 // always the same element // TODO: this should be set if the atom cannot be changed because of KEEP_NEIGHBORS
		, NO_ADDITION = 1<<2 // no atoms can be added
		, NO_REMOVAL = 1<<3 // atom cannot be removed // TODO: this should be set if the atom cannot be removed because of KEEP_NEIGHBORS
		, KEEP_NEIGHBORS = 1<<4 // the neighboring atoms must remain (bond type can change) // TODO: this should be automatically set if KEEP_NEIGHBORS_AND_BONDS is specified
		, KEEP_NEIGHBORS_AND_BONDS = 1<<5 // the neighboring atoms must remain (bond type cannot change)
		//, KEEP_BONDS = 1 << 6 // keep existing bonds as they are (neighbors can change, but not the bond order) TODO: implement
		, FULL_LOCK = NO_MUTATION | NO_ADDITION | NO_REMOVAL | KEEP_NEIGHBORS | KEEP_NEIGHBORS_AND_BONDS
	};

	static const std::vector<LockingMask> atom_locks;
	static std::string lockToString(int lock);
	static std::vector<std::string> lockingMaskToString(int mask);

	MolpherAtom(RDKit::Atom* atom);
	MolpherAtom(std::string symbol);
	MolpherAtom(std::string symbol, int formal_charge);
	MolpherAtom(const MolpherAtom& other);
	~MolpherAtom();

	MolpherAtom& operator=(const MolpherAtom&);

	bool isLocked() const;
	RDKit::Atom* asRDAtom() const;

	int getLockingMask() const;
	unsigned int getAtomicNum() const;
	double getMass() const;
	int getFormalCharge() const;
	std::string getSymbol() const;

	void setLockingMask(int mask);
	void setFormalCharge(int charge);
};

#endif //MOLPHER_LIB_MOLPHERATOM_HPP

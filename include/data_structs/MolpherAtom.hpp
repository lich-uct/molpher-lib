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
		UNLOCKED = 0
		, NO_MUTATION = 1<<1
		, NO_ADDITION = 1<<2
		, NO_REMOVAL = 1<<3
		, KEEP_NEIGHBORS = 1<<4
		, FULL_LOCK = NO_MUTATION | NO_ADDITION | NO_REMOVAL | KEEP_NEIGHBORS
	};

	static const std::vector<LockingMask> atom_locks;
	static std::string lockToString(int lock);

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

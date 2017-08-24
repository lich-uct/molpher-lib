//
// Created by sichom on 7/26/17.
//

#include <core/misc/SynchRand.h>
#include "AtomLibraryImpl.hpp"

std::unique_ptr<AtomLibrary> AtomLibrary::AtomLibraryImpl::default_lib(
		new AtomLibrary(std::vector<std::shared_ptr<MolpherAtom>>(
		{
				std::make_shared<MolpherAtom>(MolpherAtom("C"))
				, std::make_shared<MolpherAtom>(MolpherAtom("O"))
				, std::make_shared<MolpherAtom>(MolpherAtom("S"))
				, std::make_shared<MolpherAtom>(MolpherAtom("N"))
				, std::make_shared<MolpherAtom>(MolpherAtom("F"))
				, std::make_shared<MolpherAtom>(MolpherAtom("Cl"))
				, std::make_shared<MolpherAtom>(MolpherAtom("Br"))
				, std::make_shared<MolpherAtom>(MolpherAtom("I"))
		}))
);

AtomLibrary::AtomLibrary(const std::vector<std::shared_ptr<MolpherAtom>>& atoms) :
pimpl(new AtomLibraryImpl(atoms))
{
	// no action
}

const AtomLibrary &AtomLibrary::getDefaultLibrary() {
	return AtomLibraryImpl::getDefaultLibrary();
}

void AtomLibrary::setDefaultLibrary(const AtomLibrary &new_default) {
	AtomLibraryImpl::setDefaultLibrary(new_default);
}

const MolpherAtom &AtomLibrary::getRandomAtom() const {
	return pimpl->getRandomAtom();
}

AtomLibrary::AtomLibrary(const AtomLibrary &other)
: pimpl(new AtomLibraryImpl(other.getAtoms()))
{
	// no action
}

AtomLibrary::~AtomLibrary() = default;

AtomLibrary &AtomLibrary::operator=(const AtomLibrary &other) {
	pimpl = std::move(std::unique_ptr<AtomLibraryImpl>(new AtomLibraryImpl(other.getAtoms())));
	return *this;
}

std::vector<std::shared_ptr<MolpherAtom>> AtomLibrary::getAtoms() const {
	return pimpl->getAtoms();
}

AtomLibrary::AtomLibraryImpl::AtomLibraryImpl(const std::vector<std::shared_ptr<MolpherAtom>>& atoms) {
	// FIXME: check for an empty atom list and raise exception
	this->atoms.clear();
	for (auto atom : atoms) {
		this->atoms.push_back(std::make_shared<MolpherAtom>(*atom));
	}
}

const AtomLibrary &AtomLibrary::AtomLibraryImpl::getDefaultLibrary() {
	return *default_lib;
}

void AtomLibrary::AtomLibraryImpl::setDefaultLibrary(const AtomLibrary &new_default) {
	*default_lib = new_default;
}

const MolpherAtom &AtomLibrary::AtomLibraryImpl::getRandomAtom() const {
	int idx = SynchRand::GetRandomNumber( (unsigned int) atoms.size() - 1);
	return *atoms[idx];
}

std::vector<std::shared_ptr<MolpherAtom>> AtomLibrary::AtomLibraryImpl::getAtoms() const {
	std::vector<std::shared_ptr< MolpherAtom>> ret;
	for (auto atom : atoms) {
		ret.push_back(std::make_shared<MolpherAtom>(*atom));
	}
	return ret;
}
//
// Created by sichom on 7/26/17.
//

#include <core/misc/SynchRand.h>

#include "AtomLibraryImpl.hpp"

std::unique_ptr<AtomLibrary> AtomLibrary::AtomLibraryImpl::default_lib(
        new AtomLibrary(
        		std::vector<std::shared_ptr<MolpherAtom>>(
						{
								std::make_shared<MolpherAtom>(MolpherAtom("C"))
								, std::make_shared<MolpherAtom>(MolpherAtom("O"))
								, std::make_shared<MolpherAtom>(MolpherAtom("N"))
								, std::make_shared<MolpherAtom>(MolpherAtom("F"))
								, std::make_shared<MolpherAtom>(MolpherAtom("S"))
								, std::make_shared<MolpherAtom>(MolpherAtom("Cl"))
								, std::make_shared<MolpherAtom>(MolpherAtom("Br"))
								, std::make_shared<MolpherAtom>(MolpherAtom("I"))
						}),
				std::vector<double>(
						{
								73.12,
								11.741,
								11.318,
								1.379,
								1.295,
								0.8195,
								0.1703,
								0.01632
						})
        )
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

AtomLibrary::AtomLibrary(const std::vector<std::shared_ptr<MolpherAtom>>& atoms, const std::vector<double>& atom_probabilities)
        : pimpl(new AtomLibraryImpl(atoms, atom_probabilities))
{
    // no action
}

AtomLibrary::~AtomLibrary() = default;

AtomLibrary &AtomLibrary::operator=(const AtomLibrary &other) {
	pimpl = std::move(std::unique_ptr<AtomLibraryImpl>(new AtomLibraryImpl(other.getAtoms(), other.getAtomProbabilities())));
	return *this;
}

std::vector<std::shared_ptr<MolpherAtom>> AtomLibrary::getAtoms() const {
	return pimpl->getAtoms();
}

AtomLibrary::AtomLibrary(const MolpherMol &mol) : pimpl(new AtomLibraryImpl(mol)) {
	// no action
}

std::vector<double> AtomLibrary::getAtomProbabilities() const {
	return pimpl->getAtomProbabilities();
}

void AtomLibrary::setAtomProbabilities(const std::vector<double> &atom_probabilities) {
	return pimpl->setAtomProbabilities(atom_probabilities);
}

// impl stuff:

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
    if (atom_probabilities.empty()){
        int idx = SynchRand::GetRandomNumber( (unsigned int) atoms.size() - 1);
        return *atoms[idx];
    }else {
        double sum = 0;
        for (int i = 0; i < atom_probabilities.size(); i++) {
            sum += atom_probabilities[i];
        }
        double rnd = (double) SynchRand::GetRandomNumber(RAND_MAX) / RAND_MAX * sum;
        for (int idx = 0; idx < atom_probabilities.size(); idx++) {
            if (rnd < atom_probabilities[idx]) {
                return *atoms[idx];
            }
            rnd -= atom_probabilities[idx];
        }
    }
}

std::vector<std::shared_ptr<MolpherAtom>> AtomLibrary::AtomLibraryImpl::getAtoms() const {
	std::vector<std::shared_ptr< MolpherAtom>> ret;
	for (auto atom : atoms) {
		ret.push_back(std::make_shared<MolpherAtom>(*atom));
	}
	return ret;
}

AtomLibrary::AtomLibraryImpl::AtomLibraryImpl(const MolpherMol &mol)
: AtomLibraryImpl(mol.getAtoms())
{
	// no action
}

AtomLibrary::AtomLibraryImpl::AtomLibraryImpl(const std::vector<std::shared_ptr<MolpherAtom>> &atoms,
                                              const std::vector<double> &atom_probabilities)
        : AtomLibraryImpl(atoms)
{
    this->atom_probabilities.clear();
    for (auto prob : atom_probabilities) {
        this->atom_probabilities.push_back(prob);
    }
}

std::vector<double> AtomLibrary::AtomLibraryImpl::getAtomProbabilities() const {
    std::vector<double> probabilities;
    for (auto prob : atom_probabilities){
        probabilities.push_back(prob);
    }
    return probabilities;
}

void
AtomLibrary::AtomLibraryImpl::setAtomProbabilities(const std::vector<double> &atom_probabilities) {
	this->atom_probabilities = atom_probabilities;
}


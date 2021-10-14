//
// Created by sichom on 7/26/17.
//

#ifndef MOLPHER_LIB_ATOMLIBRARYIMPL_HPP
#define MOLPHER_LIB_ATOMLIBRARYIMPL_HPP

#include "morphing/AtomLibrary.hpp"


class AtomLibrary::AtomLibraryImpl {

private:
	std::vector<std::shared_ptr<const MolpherAtom>> atoms;
    std::vector<double> atom_probabilities;
	static std::unique_ptr<AtomLibrary> default_lib;

public:
	AtomLibraryImpl(const std::vector<std::shared_ptr<MolpherAtom>>& atoms);

    AtomLibraryImpl(const std::vector<std::shared_ptr<MolpherAtom>>& atoms, const std::vector<double>& atom_probabilities);
    AtomLibraryImpl(const MolpherMol& mol);

	static const AtomLibrary& getDefaultLibrary();
	static void setDefaultLibrary(const AtomLibrary& new_default);

    const MolpherAtom& getRandomAtom() const;
	std::vector<std::shared_ptr<MolpherAtom>> getAtoms() const;
    std::vector<double> getAtomProbabilities() const;
	void setAtomProbabilities(const std::vector<double>& atom_probabilities);
};

#endif //MOLPHER_LIB_ATOMLIBRARYIMPL_HPP

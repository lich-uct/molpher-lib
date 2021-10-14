//
// Created by sichom on 7/26/17.
//

#ifndef MOLPHER_LIB_ATOMLIBRARY_HPP
#define MOLPHER_LIB_ATOMLIBRARY_HPP

#include "data_structs/MolpherMol.hpp"
#include <vector>


class AtomLibrary {

private:
	class AtomLibraryImpl;
	std::unique_ptr<AtomLibraryImpl> pimpl;

public:
	AtomLibrary(const std::vector<std::shared_ptr<MolpherAtom>>& atoms);
    AtomLibrary(const std::vector<std::shared_ptr<MolpherAtom>>& atoms, const std::vector<double>& atom_probabilities);
    AtomLibrary(const MolpherMol& mol);
	AtomLibrary(const AtomLibrary& other);
	~AtomLibrary();

	AtomLibrary& operator=(const AtomLibrary&);

	static const AtomLibrary& getDefaultLibrary();
	static void setDefaultLibrary(const AtomLibrary& new_default);

    const MolpherAtom& getRandomAtom() const;
	std::vector<std::shared_ptr<MolpherAtom>> getAtoms() const;
    std::vector<double> getAtomProbabilities() const;
	void setAtomProbabilities(const std::vector<double>& atom_probabilities);

};

#endif //MOLPHER_LIB_ATOMLIBRARY_HPP

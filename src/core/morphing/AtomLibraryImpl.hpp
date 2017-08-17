//
// Created by sichom on 7/26/17.
//

#ifndef MOLPHER_LIB_ATOMLIBRARYIMPL_HPP
#define MOLPHER_LIB_ATOMLIBRARYIMPL_HPP

#include "morphing/AtomLibrary.hpp"

class AtomLibrary::AtomLibraryImpl {

private:
	std::vector<std::shared_ptr<const MolpherAtom>> atoms;
	static std::unique_ptr<AtomLibrary> default_lib;

public:
	AtomLibraryImpl(const std::vector<std::shared_ptr<MolpherAtom>>& atoms);

	static const AtomLibrary& getDefaultLibrary();
	static void setDefaultLibrary(const AtomLibrary& new_default);

	const MolpherAtom& getRandomAtom() const;
	std::vector<std::shared_ptr<MolpherAtom>> getAtoms() const;

};

#endif //MOLPHER_LIB_ATOMLIBRARYIMPL_HPP

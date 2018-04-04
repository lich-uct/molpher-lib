//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_ADDATOM_HPP
#define MOLPHER_LIB_ADDATOM_HPP

#include <morphing/AtomLibrary.hpp>
#include "morphing/operators/MorphingOperator.hpp"

class AddAtom : public MorphingOperator {
public:
	AddAtom();
	AddAtom(const AtomLibrary& atom_library);

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();
	virtual std::string getName() const;

	const std::vector<unsigned int>& getOpenIndices();
	std::vector<std::shared_ptr<MolpherAtom>> getOpenAtoms();

private:
	class AddAtomImpl;
	std::shared_ptr<AddAtomImpl> pimpl;
};

#endif //MOLPHER_LIB_ADDATOM_HPP

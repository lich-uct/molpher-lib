//
// Created by sichom on 8/28/17.
//

#ifndef MOLPHER_LIB_REMOVEATOM_HPP
#define MOLPHER_LIB_REMOVEATOM_HPP


#include "MorphingOperator.hpp"

class RemoveAtom : public MorphingOperator {
public:
	RemoveAtom();

	virtual void setOriginal(std::shared_ptr<MolpherMol> mol);
	virtual std::shared_ptr<MolpherMol> morph();
	virtual std::string getName() const;

	const std::vector<unsigned int>& getMarkedIndices();
	std::vector<std::shared_ptr<MolpherAtom>> getMarkedAtoms();

private:
	class RemoveAtomImpl;
	std::shared_ptr<RemoveAtomImpl> pimpl;

};


#endif //MOLPHER_LIB_REMOVEATOM_HPP

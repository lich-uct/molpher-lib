//
// Created by sichom on 7/25/17.
//

#ifndef MOLPHER_LIB_MORPHINGOPERATOR_HPP
#define MOLPHER_LIB_MORPHINGOPERATOR_HPP

#include "data_structs/MolpherMol.hpp"

class MorphingOperator {
public:
	MorphingOperator();
	virtual ~MorphingOperator();

	std::shared_ptr<MolpherMol> getOriginal();
	virtual void setOriginal(std::shared_ptr<MolpherMol>) = 0;
	virtual std::shared_ptr<MolpherMol> morph() = 0; // FIXME: make this a const method to guarantee thread safety

protected:
	class MorphingOperatorImpl;
	void setMorphingOperatorPimpl(std::shared_ptr<MorphingOperatorImpl> pimpl);

private:
	std::shared_ptr<MorphingOperatorImpl> pimpl;

};

#endif //MOLPHER_LIB_MORPHINGOPERATOR_HPP

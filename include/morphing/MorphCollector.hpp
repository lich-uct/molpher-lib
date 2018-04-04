//
// Created by sichom on 3/26/18.
//

#ifndef MOLPHER_LIB_MORPHCOLLECTOR_HPP
#define MOLPHER_LIB_MORPHCOLLECTOR_HPP

#include <memory>
#include <morphing/operators/MorphingOperator.hpp>
#include "data_structs/MolpherMol.hpp"

class MorphCollector {

private:
	class MorphCollectorImpl;
	std::unique_ptr<MorphCollectorImpl> pimpl;

public:
	MorphCollector();
	virtual ~MorphCollector();
	virtual void operator()(std::shared_ptr<MolpherMol> morph, std::shared_ptr<MorphingOperator> operator_) = 0;

};

#endif //MOLPHER_LIB_MORPHCOLLECTOR_HPP

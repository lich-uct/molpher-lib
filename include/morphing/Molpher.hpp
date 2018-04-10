//
// Created by sichom on 8/30/17.
//

#ifndef MOLPHER_LIB_MOLPHER_HPP
#define MOLPHER_LIB_MOLPHER_HPP

#include <memory>
#include <vector>
#include "morphing/operators/MorphingOperator.hpp"
#include "data_structs/MolpherMol.hpp"
#include "MorphCollector.hpp"

class Molpher {
private:
	class MolpherImpl;
	std::unique_ptr<MolpherImpl> pimpl;

public:
	Molpher(
			std::shared_ptr<MolpherMol> mol
			, const std::vector<std::shared_ptr<MorphingOperator> >& operators
			, unsigned int threads
			, unsigned int attempts
	);
	Molpher(
			std::shared_ptr<MolpherMol> mol
			, const std::vector<std::shared_ptr<MorphingOperator> >& operators
			, unsigned int threads
			, unsigned int attempts
			, const std::vector<std::shared_ptr<MorphCollector> >& collectors
	);
	~Molpher();

	void operator()();
	std::vector<std::shared_ptr<MolpherMol> > getMorphs();
	std::shared_ptr<MolpherMol> getOriginal();
	void reset(std::shared_ptr<MolpherMol> original);
};

#endif //MOLPHER_LIB_MOLPHER_HPP

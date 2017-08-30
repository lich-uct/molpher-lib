//
// Created by sichom on 8/30/17.
//

#ifndef MOLPHER_LIB_MOLPHERIMPL_HPP
#define MOLPHER_LIB_MOLPHERIMPL_HPP

#include "core/misc/global_types.h"
#include "morphing/Molpher.hpp"
#include "MorphCalculator.hpp"

class Molpher::MolpherImpl {
private:
	std::vector<MorphingOperator*> operators;
	MorphCalculator calc;
	std::shared_ptr<MolpherMol> original;
	unsigned int threads;
	unsigned int attempts;

	ConcurrentMolVector morphs;
	tbb::atomic<unsigned int> failures;
	tbb::atomic<unsigned int> empty_mols;
public:
	MolpherImpl(std::shared_ptr<MolpherMol> mol, const std::vector<MorphingOperator*>& operators, unsigned int threads, unsigned int attempts);

	void operator()();
	std::vector<std::shared_ptr<MolpherMol> > getMorphs();
	std::shared_ptr<MolpherMol> getOriginal();
	void reset(std::shared_ptr<MolpherMol> original);
};

#endif //MOLPHER_LIB_MOLPHERIMPL_HPP

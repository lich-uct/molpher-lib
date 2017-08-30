//
// Created by sichom on 8/30/17.
//

#ifndef MOLPHER_LIB_MORPHCALCULATOR_HPP
#define MOLPHER_LIB_MORPHCALCULATOR_HPP


#include <boost/shared_ptr.hpp>
#include <tbb/atomic.h>
#include <tbb/blocked_range.h>
#include "morphing/operators/MorphingOperator.hpp"

class MorphCalculator {
public:
	MorphCalculator(
			std::vector<MorphingOperator*>& operators,
			ConcurrentMolVector& morphs,
			tbb::atomic<unsigned int>& failures,
			tbb::atomic<unsigned int>& empty_mols
	);

	void operator()(const tbb::blocked_range<int> &r) const;

private:
	std::vector<MorphingOperator*>& operators;
	ConcurrentMolVector& morphs;
	tbb::atomic<unsigned int>& mMorphingFailureCount;
	tbb::atomic<unsigned int>& mMorphEmptyCount;
};


#endif //MOLPHER_LIB_MORPHCALCULATOR_HPP

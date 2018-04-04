//
// Created by sichom on 8/30/17.
//

#ifndef MOLPHER_LIB_MORPHCALCULATOR_HPP
#define MOLPHER_LIB_MORPHCALCULATOR_HPP


#include <boost/shared_ptr.hpp>
#include <tbb/atomic.h>
#include <tbb/blocked_range.h>
#include <morphing/MorphCollector.hpp>
#include "morphing/operators/MorphingOperator.hpp"

class MorphCalculator {
public:
	MorphCalculator(
			std::vector<std::shared_ptr<MorphingOperator> >& operators,
			ConcurrentMolVector& morphs,
			tbb::atomic<unsigned int>& failures,
			tbb::atomic<unsigned int>& empty_mols
	);
	MorphCalculator(
			std::vector<std::shared_ptr<MorphingOperator> >& operators,
			ConcurrentMolVector& morphs,
			tbb::atomic<unsigned int>& failures,
			tbb::atomic<unsigned int>& empty_mols,
			std::shared_ptr<MorphCollector> collector
	);

	void operator()(const tbb::blocked_range<int> &r) const;

private:
	std::vector<std::shared_ptr<MorphingOperator> >& operators;
	ConcurrentMolVector& morphs;
	tbb::atomic<unsigned int>& mMorphingFailureCount;
	tbb::atomic<unsigned int>& mMorphEmptyCount;
	std::shared_ptr<MorphCollector> collector;
};


#endif //MOLPHER_LIB_MORPHCALCULATOR_HPP

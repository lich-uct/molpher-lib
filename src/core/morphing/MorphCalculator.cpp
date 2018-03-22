//
// Created by sichom on 8/30/17.
//

#include <core/chem/ChemicalAuxiliary.h>
#include "core/misc/SynchRand.h"
#include "MorphCalculator.hpp"
#include "core/misc/inout.h"

MorphCalculator::MorphCalculator(std::vector<std::shared_ptr<MorphingOperator> > &operators,
								 ConcurrentMolVector& morphs,
								 tbb::atomic<unsigned int>& failures,
								 tbb::atomic<unsigned int>& empty_mols
)
:
operators(operators)
, morphs(morphs)
, mMorphEmptyCount(failures)
, mMorphingFailureCount(empty_mols)
{
	// no action
}

void MorphCalculator::operator()(const tbb::blocked_range<int> &r) const {
	for (int i = r.begin(); i != r.end(); ++i) {
		int randPos = SynchRand::GetRandomNumber(operators.size() - 1);
		std::shared_ptr<MorphingOperator> operator_ = operators[randPos];

		std::shared_ptr<MolpherMol> new_mol(nullptr);
		try {
			new_mol = operator_->morph();
		} catch (const std::exception &exc) {
			SynchCerr("Morphing failure due to an error: " + std::string(exc.what()));
			++mMorphingFailureCount; // atomic
			return;
		}

		if (new_mol) {
			morphs.push_back(new_mol);
		} else {
			SynchCerr("Morphing failure: morphing method returned empty molecule");
			++mMorphEmptyCount;
		}
	}
}

//
// Created by sichom on 8/30/17.
//

#include <core/chem/ChemicalAuxiliary.h>
#include "core/misc/SynchRand.h"
#include "MorphCalculator.hpp"
#include "core/misc/inout.h"

MorphCalculator::MorphCalculator(
		std::vector<std::shared_ptr<MorphingOperator> > &operators
		, ConcurrentMolVector &morphs
		, tbb::atomic<unsigned int> &failures
		, tbb::atomic<unsigned int> &empty_mols
		, const std::vector<std::shared_ptr<MorphCollector> >& collectors
) :
operators(operators)
, morphs(morphs)
, mMorphingFailureCount(failures)
, mMorphEmptyCount(empty_mols)
, collectors(collectors)
{
	// no action
}

MorphCalculator::MorphCalculator(
		std::vector<std::shared_ptr<MorphingOperator> > &operators,
		ConcurrentMolVector& morphs,
		tbb::atomic<unsigned int>& failures,
		tbb::atomic<unsigned int>& empty_mols
)
:
MorphCalculator::MorphCalculator(
		operators
		, morphs
		, failures
		, empty_mols
		, std::vector<std::shared_ptr<MorphCollector>>())
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
			SynchCerr(
					"Morphing failure: " + std::string(exc.what())
					+ "\n\tParent molecule: " + operator_->getOriginal()->getSMILES()
					+ "\n\tParent operator: " + operator_->getName()
			);
			++mMorphingFailureCount; // atomic
			continue;
		}

		if (!new_mol) {
//			SynchCerr("Morphing method returned empty molecule.");
			++mMorphEmptyCount;
			continue;
		}

		try {
			for (auto collector : collectors) {
				(*collector)(new_mol, operator_);
			}
		} catch (const std::exception &exc) {
			SynchCerr(
					"Failed to collect morph due to: " + std::string(exc.what())
					+ "\n\tParent molecule: " + operator_->getOriginal()->getSMILES()
					+ "\n\tParent operator: " + operator_->getName()
					+ "\n\tGenerated morph: " + new_mol->getSMILES()
			);
			continue;
		}

		morphs.push_back(new_mol);
	}
}

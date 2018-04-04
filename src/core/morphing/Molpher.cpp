//
// Created by sichom on 8/30/17.
//

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <core/misc/inout.h>
#include "MolpherImpl.hpp"

Molpher::Molpher(std::shared_ptr<MolpherMol> mol, const std::vector<std::shared_ptr<MorphingOperator> > &operators, unsigned int threads,
				 unsigned int attempts)
:
pimpl(new MolpherImpl(mol, operators, threads, attempts))
{
	// no action
}

void Molpher::operator()() {
	(*pimpl)();
}

std::vector<std::shared_ptr<MolpherMol> > Molpher::getMorphs() {
	return pimpl->getMorphs();
}

std::shared_ptr<MolpherMol> Molpher::getOriginal() {
	return pimpl->getOriginal();
}

void Molpher::reset(std::shared_ptr<MolpherMol> original) {
	pimpl->reset(original);
}

Molpher::Molpher(std::shared_ptr<MolpherMol> mol, const std::vector<std::shared_ptr<MorphingOperator> > &operators,
				 unsigned int threads, unsigned int attempts, std::shared_ptr<MorphCollector> collector)
: pimpl(new MolpherImpl(mol, operators, threads, attempts, collector))
{
	// no action
}

Molpher::~Molpher() = default;

// implementation

Molpher::MolpherImpl::MolpherImpl(std::shared_ptr<MolpherMol> mol,
								  const std::vector<std::shared_ptr<MorphingOperator> > &operators,
								  unsigned int threads, unsigned int attempts, std::shared_ptr<MorphCollector> collector)
		:
		original(mol)
		, operators(operators)
		, threads(threads)
		, attempts(attempts)
		, morphs()
		, failures(0)
		, empty_mols(0)
		, calc(this->operators, morphs, failures, empty_mols, collector)
{
	for (const auto& oper : operators) {
		try {
			oper->setOriginal(original);
		} catch (const std::exception& exc) {
			SynchCerr(
					"Failed to set original. The following error occured: " + std::string(exc.what())
					+ "\n\tOriginal: " + original->getSMILES()
					+ "\n\tFailed operator: " + oper->getName()
					+ "\n\tParent: " + original->getParentSMILES()
					+ "\n\tParent operator: " + original->getParentOper()
					, "ERROR: "
			);
			throw exc;
		}
	}
}

Molpher::MolpherImpl::MolpherImpl(std::shared_ptr<MolpherMol> mol, const std::vector<std::shared_ptr<MorphingOperator> > &operators,
								  unsigned int threads, unsigned int attempts)
:
MolpherImpl(mol, operators, threads, attempts, nullptr)
{
	// no action
}

void Molpher::MolpherImpl::operator()() {
	tbb::task_group_context tbbCtx;
	tbb::task_scheduler_init scheduler;
	if (threads > 0) {
		scheduler.terminate();
		scheduler.initialize(threads);
	}

	try {
		tbb::parallel_for(tbb::blocked_range<int>(0, attempts),
					  calc, tbb::auto_partitioner(), tbbCtx);
	} catch (const std::exception& exp) {
		Cerr(
			"Morphing finished prematurely due to an exception: " + std::string(exp.what())
			+ "\n\tOriginal molecule: " + original->getSMILES()
		);
	}
}

std::vector<std::shared_ptr<MolpherMol> > Molpher::MolpherImpl::getMorphs() {
	std::vector<std::shared_ptr<MolpherMol> > ret;
	for (auto mol : morphs) {
		ret.push_back(mol);
	}
	morphs.clear();
	return ret;
}

std::shared_ptr<MolpherMol> Molpher::MolpherImpl::getOriginal() {
	return original;
}

void Molpher::MolpherImpl::reset(std::shared_ptr<MolpherMol> original) {
	this->original = original;
	for (auto oper : operators) {
		oper->setOriginal(this->original);
	}
	morphs.clear();
	failures = 0;
	empty_mols = 0;
}

//
// Created by sichom on 8/30/17.
//

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
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

Molpher::~Molpher() = default;

// implementation

Molpher::MolpherImpl::MolpherImpl(std::shared_ptr<MolpherMol> mol, const std::vector<std::shared_ptr<MorphingOperator> > &operators,
								  unsigned int threads, unsigned int attempts)
:
original(mol)
, operators(operators)
, threads(threads)
, attempts(attempts)
, morphs()
, failures(0)
, empty_mols(0)
, calc(this->operators, morphs, failures, empty_mols)
{
	for (auto oper : operators) {
		oper->setOriginal(original);
	}
}

void Molpher::MolpherImpl::operator()() {
	tbb::task_group_context tbbCtx;
	tbb::task_scheduler_init scheduler;
	if (threads > 0) {
		scheduler.terminate();
		scheduler.initialize(threads);
	}

	tbb::parallel_for(tbb::blocked_range<int>(0, attempts),
					  calc, tbb::auto_partitioner(), tbbCtx);
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

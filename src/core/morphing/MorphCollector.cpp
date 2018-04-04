//
// Created by sichom on 3/26/18.
//

#include "MorphCollectorImpl.hpp"

MorphCollector::MorphCollector()
: pimpl(new MorphCollectorImpl())
{
	// no action
}

MorphCollector::~MorphCollector() = default;

MorphCollector::MorphCollectorImpl::MorphCollectorImpl() = default;

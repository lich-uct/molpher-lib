/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "core/misc/inout.h"
#include "core/misc/SynchRand.h"

#include "operations/FilterMorphsOper.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "FilterMorphsOperImpl.hpp"

FilterMorphsOper::FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, bool verbose) : 
pimpl(new FilterMorphsOper::FilterMorphsOperImpl(expTree, verbose))
{
    setTreeOperPimpl(pimpl);
}

FilterMorphsOper::FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, FilterMorphsOper::MorphFilters filters, bool verbose) : 
pimpl(new FilterMorphsOper::FilterMorphsOperImpl(expTree, static_cast<int>(filters), verbose))
{
    setTreeOperPimpl(pimpl);
}

FilterMorphsOper::FilterMorphsOper(bool verbose) : 
pimpl(new FilterMorphsOper::FilterMorphsOperImpl(verbose))
{
    setTreeOperPimpl(pimpl);
}

FilterMorphsOper::FilterMorphsOper(FilterMorphsOper::MorphFilters filters, bool verbose) : 
pimpl(new FilterMorphsOper::FilterMorphsOperImpl(static_cast<int>(filters), verbose))
{
    setTreeOperPimpl(pimpl);
}

void FilterMorphsOper::operator()() {
    (*pimpl)();
}

// pimpl

FilterMorphsOper::FilterMorphsOperImpl::FilterMorphsOperImpl(
std::shared_ptr<ExplorationTree> expTree
, bool verbose
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, filters(FilterMorphsOper::MorphFilters::ALL)
, verbose(verbose)
{
    // no action
}

FilterMorphsOper::FilterMorphsOperImpl::FilterMorphsOperImpl(
int filters
, bool verbose
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, filters(filters | FilterMorphsOper::MorphFilters::DUPLICATES)
, verbose(verbose)
{
    // no action
}

FilterMorphsOper::FilterMorphsOperImpl::FilterMorphsOperImpl(bool verbose) : 
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, filters(FilterMorphsOper::MorphFilters::ALL)
, verbose(verbose)
{
    // no action
}

FilterMorphsOper::FilterMorphsOperImpl::FilterMorphsOperImpl(
std::shared_ptr<ExplorationTree> expTree
, int filters
, bool verbose
) 
:
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, filters(filters | FilterMorphsOper::MorphFilters::DUPLICATES)
, verbose(verbose)
{
    // no action
}

FilterMorphsOper::FilterMorphsOperImpl::FilterMorphs::FilterMorphs(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_impl,
        size_t globalMorphCount, ConcurrentMolVector &morphs, std::vector<bool> &survivors, int filters, bool verbose
        ) :
mTreePimpl(tree_impl),
mGlobalMorphCount(globalMorphCount),
mMorphs(morphs),
mSurvivors(survivors),
mVerboseOutput(verbose),
mFilters(filters){
    assert(mMorphs.size() == mSurvivors.size());
}

void FilterMorphsOper::FilterMorphsOperImpl::FilterMorphs::operator()(const tbb::blocked_range<size_t> &r) const {

    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        
        std::string morph_smiles = mMorphs[idx]->getSMILES();

        bool mightSurvive = true;
        if (mFilters & MorphFilters::PROBABILITY) {
            double acceptProbability = 1.0;
            bool isTarget = (morph_smiles == mTreePimpl->target.getSMILES());
            if (idx >= mTreePimpl->params.cntCandidatesToKeep && !isTarget) {
                acceptProbability =
                        0.25 - (idx - mTreePimpl->params.cntCandidatesToKeep) /
                        ((mGlobalMorphCount - mTreePimpl->params.cntCandidatesToKeep) * 4.0);
            }
            mightSurvive = SynchRand::GetRandomNumber(0, 99) < (int) (acceptProbability * 100);
        }
        
        if (mightSurvive) {
            bool isDead = false;
            bool badWeight = false;
            bool badSascore = false;
            bool alreadyInTree = false;
            bool alreadyTriedByParent = false;
            bool tooManyProducedMorphs = false;

            // Tests are ordered according to their cost.

            if (mFilters & MorphFilters::WEIGHT) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    double weight = mMorphs[idx]->getMolecularWeight();
                    badWeight =
                            ( weight <
                            mTreePimpl->params.minAcceptableMolecularWeight) ||
                            (weight >
                            mTreePimpl->params.maxAcceptableMolecularWeight);
                    if (badWeight && mVerboseOutput) {
                        std::stringstream ss;
                        ss << "bad weight: " << morph_smiles << " : " << weight;
                        SynchCout(ss.str());
                    }
                }
            }

            if (mFilters & MorphFilters::DUPLICATES) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    TreeMap::const_accessor ac;
                    if (mTreePimpl->treeMap.find(ac, morph_smiles)) {
                        alreadyInTree = true;
                        SynchCout("Duplicate molecule found: " + morph_smiles);
                    }
                }
            }

            if (mFilters & MorphFilters::HISTORIC_DESCENDENTS) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    TreeMap::const_accessor ac;
                    if (mTreePimpl->treeMap.find(ac, mMorphs[idx]->getParentSMILES())) {
                        auto descendants = ac->second->getHistoricDescendants();
                        alreadyTriedByParent = (
                                descendants.find(morph_smiles)
                                !=
                                descendants.end());
                    } else {
                        assert(false);
                    }
                }
            }

            if (mFilters & MorphFilters::MAX_DERIVATIONS) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    MorphDerivationMap::const_accessor ac;
                    if (mTreePimpl->morphDerivations.find(ac, morph_smiles)) {
                        tooManyProducedMorphs =
                                (ac->second > mTreePimpl->params.cntMaxMorphs);
                    }
                    if (tooManyProducedMorphs && mVerboseOutput) {
                        std::stringstream ss;
                        ss << "too many morphs: " << morph_smiles << " : " << ac->second;
                        SynchCout(ss.str());
                    }
                }
            }
            
            if (mFilters & MorphFilters::SYNTHESIS) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    double sascore = mMorphs[idx]->getSAScore();
                    badSascore = sascore > 6.0; // questionable, it is recommended value from Ertl
                    // in case of badSascore print message
                    if (badSascore && mVerboseOutput) {
                        std::stringstream ss;
                        ss << "bad SAScore: " << morph_smiles << " : " << sascore;
                        SynchCout(ss.str());
                    }
                }
            }

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs);
            mSurvivors[idx] = !isDead;
        } else {
            mSurvivors[idx] = false;
            if (mVerboseOutput) {
                std::stringstream ss;
                ss << "probability filtered: " << morph_smiles;
                SynchCout(ss.str());
            }
        }
    }
}

void FilterMorphsOper::FilterMorphsOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_impl = tree->pimpl;
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (tree_impl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree_impl->threadCnt);
        }

        assert(tree_impl->candidates.size() == tree_impl->candidatesMask.size());
        FilterMorphs filterMorphs(tree_impl, tree_impl->candidates.size(), tree_impl->candidates, tree_impl->candidatesMask, filters, verbose);
        tbb::parallel_for(
                tbb::blocked_range<size_t>(0, tree_impl->candidates.size()),
                filterMorphs, tbb::auto_partitioner(), tbbCtx);
    } else {
        throw std::runtime_error("Cannot filter morphs. No tree associated with this instance.");
    }
}



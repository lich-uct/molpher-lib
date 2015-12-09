
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"

#include "molpher_API/operations/FilterMorphsOper.hpp"

FilterMorphsOper::FilterMorphsOper(ExplorationTree& expTree) : TreeOperation(expTree), filters(MorphFilters::ALL) {
    // no action
}

FilterMorphsOper::FilterMorphsOper(ExplorationTree& expTree, int filters) : TreeOperation(expTree), filters(filters | MorphFilters::DUPLICATES) {
    // no action
}

FilterMorphsOper::FilterMorphs::FilterMorphs(PathFinderContext &ctx,
        size_t globalMorphCount, ExplorationTree::MoleculeVector &morphs, std::vector<bool> &survivors, int filters
        ) :
mCtx(ctx),
mGlobalMorphCount(globalMorphCount),
mMorphs(morphs),
mSurvivors(survivors),
mFilters(filters){
    assert(mMorphs.size() == mSurvivors.size());
}

void FilterMorphsOper::FilterMorphs::operator()(const tbb::blocked_range<size_t> &r) const {

    for (size_t idx = r.begin(); idx != r.end(); ++idx) {

        bool mightSurvive = true;
        if (mFilters & MorphFilters::PROBABILITY) {
            double acceptProbability = 1.0;
            bool isTarget = (mMorphs[idx].smile == mCtx.target.smile);
            if (idx >= mCtx.params.cntCandidatesToKeep && !isTarget) {
                acceptProbability =
                        0.25 - (idx - mCtx.params.cntCandidatesToKeep) /
                        ((mGlobalMorphCount - mCtx.params.cntCandidatesToKeep) * 4.0);
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
            // Added test for SAScore

            if (mFilters & MorphFilters::WEIGHT) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    badWeight =
                            (mMorphs[idx].molecularWeight <
                            mCtx.params.minAcceptableMolecularWeight) ||
                            (mMorphs[idx].molecularWeight >
                            mCtx.params.maxAcceptableMolecularWeight);
                    if (badWeight) {
                        std::stringstream ss;
                        ss << "bad weight: " << mMorphs[idx].smile << " : " << mMorphs[idx].molecularWeight;
                        SynchCout(ss.str());
                    }
                }
            }

            if (mFilters & MorphFilters::DUPLICATES) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    PathFinderContext::CandidateMap::const_accessor ac;
                    if (mCtx.candidates.find(ac, mMorphs[idx].smile)) {
                        alreadyInTree = true;
                    }
                }
            }

            if (mFilters & MorphFilters::HISTORIC_DESCENDANTS) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    PathFinderContext::CandidateMap::const_accessor ac;
                    if (mCtx.candidates.find(ac, mMorphs[idx].parentSmile)) {
                        alreadyTriedByParent = (
                                ac->second.historicDescendants.find(mMorphs[idx].smile)
                                !=
                                ac->second.historicDescendants.end());
                    } else {
                        assert(false);
                    }
                }
            }

            if (mFilters & MorphFilters::COUNT) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    PathFinderContext::MorphDerivationMap::const_accessor ac;
                    if (mCtx.morphDerivations.find(ac, mMorphs[idx].smile)) {
                        tooManyProducedMorphs =
                                (ac->second > mCtx.params.cntMaxMorphs);
                    }
                    if (tooManyProducedMorphs) {
                        std::stringstream ss;
                        ss << "too many morphs: " << mMorphs[idx].smile << " : " << ac->second;
                        SynchCout(ss.str());
                    }
                }
            }
            
            if (mFilters & MorphFilters::SYNTHESIS) {
                isDead = (badWeight || badSascore || alreadyInTree ||
                        alreadyTriedByParent || tooManyProducedMorphs);
                if (!isDead) {
                    badSascore = mMorphs[idx].sascore > 6.0; // questionable, it is recommended value from Ertl
                    // in case of badSascore print message
                    if (badSascore) {
                        std::stringstream ss;
                        ss << "bad SAScore: " << mMorphs[idx].smile << " : " << mMorphs[idx].sascore;
                        SynchCout(ss.str());
                    }
                }
            }

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs);
            mSurvivors[idx] = !isDead;
        } else {
            mSurvivors[idx] = false;
            std::stringstream ss;
                        ss << "probability filtered: " << mMorphs[idx].smile;
                        SynchCout(ss.str());
        }
    }
}

void FilterMorphsOper::operator()() {
    tbb::task_group_context tbbCtx;
    tbb::task_scheduler_init scheduler;
    if (threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(threadCnt);
    }

    ExplorationTree::MoleculeVector& morphs = fetchGeneratedMorphs();
    PathFinderContext& context = fetchTreeContext();
    ExplorationTree::BoolVector& survivors = fetchGeneratedMorphsMask();
    assert(morphs.size() == survivors.size());
    FilterMorphs filterMorphs(context, morphs.size(), morphs, survivors, filters);
    tbb::parallel_for(
            tbb::blocked_range<size_t>(0, morphs.size()),
            filterMorphs, tbb::auto_partitioner(), tbbCtx);
}



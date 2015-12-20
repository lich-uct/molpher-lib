
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "chem/morphing/Morphing.hpp"

#include "molpher_API/operations/GenerateMorphsOper.hpp"

GenerateMorphsOper::GenerateMorphsOper(ExplorationTree& expTree) : TreeOperation(expTree) {
    // no action
}

GenerateMorphsOper::GenerateMorphsOper() : TreeOperation() {
    // no action
}

void GenerateMorphsOper::CollectMorphs::MorphCollector(MolpherMolecule *morph, void *functor) {
    GenerateMorphsOper::CollectMorphs *collect =
            (GenerateMorphsOper::CollectMorphs *) functor;
    (*collect)(*morph);
}

GenerateMorphsOper::CollectMorphs::CollectMorphs(ExplorationTree::MoleculeVector &morphs) :
mMorphs(morphs) {
    mCollectAttemptCount = 0;
}

void GenerateMorphsOper::CollectMorphs::operator()(const MolpherMolecule &morph) {
    ++mCollectAttemptCount; // atomic
    ExplorationTree::SmileSet::const_accessor dummy;
    if (mDuplicateChecker.insert(dummy, morph.smile)) {
        mMorphs.push_back(morph);
    } else {
        // ignore duplicate
    }
}

unsigned int GenerateMorphsOper::CollectMorphs::WithdrawCollectAttemptCount() {
    unsigned int ret = mCollectAttemptCount;
    mCollectAttemptCount = 0;
    return ret;
}

void GenerateMorphsOper::operator()() {
    if (this->tree) {
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(threadCnt);
        }

        ExplorationTree::MoleculePointerVector leaves;
        this->tree->fetchLeaves(leaves);
        PathFinderContext& context = fetchTreeContext();
        ExplorationTree::MoleculeVector& morphs = fetchGeneratedMorphs();
        morphs.clear();
        CollectMorphs collectMorphs(morphs);
        for (auto it : leaves) {
            MolpherMolecule &candidate = (*it);
            unsigned int morphAttempts = context.params.cntMorphs;
            if (candidate.distToTarget < context.params.distToTargetDepthSwitch) {
                morphAttempts = context.params.cntMorphsInDepth;
            }

            morphs.reserve(morphs.size() + morphAttempts);

            GenerateMorphs(
                    candidate,
                    morphAttempts,
                    context.fingerprintSelector,
                    context.simCoeffSelector,
                    context.chemOperSelectors,
                    context.target,
                    context.decoys,
                    tbbCtx,
                    &collectMorphs,
                    GenerateMorphsOper::CollectMorphs::MorphCollector);
            PathFinderContext::MorphDerivationMap::accessor ac;

            if (context.morphDerivations.find(ac, candidate.smile)) {
                ac->second += collectMorphs.WithdrawCollectAttemptCount();
            } else {
                context.morphDerivations.insert(ac, candidate.smile);
                ac->second = collectMorphs.WithdrawCollectAttemptCount();
            }
        }
        morphs.shrink_to_fit();
        ExplorationTree::BoolVector& survivors = fetchGeneratedMorphsMask();
        survivors.resize(morphs.size(), true);
    } else {
        throw std::runtime_error("Cannot generate morphs. No tree associated with this instance.");
    }
}

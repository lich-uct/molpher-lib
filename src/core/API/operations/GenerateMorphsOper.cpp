
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "core/chem/morphing/Morphing.hpp"

#include "operations/GenerateMorphsOper.hpp"
#include "GenerateMorphsOperImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "core/misc/inout.h"

GenerateMorphsOper::GenerateMorphsOper(std::shared_ptr<ExplorationTree> expTree, bool set_tree_ownership) : 
pimpl(new GenerateMorphsOper::GenerateMorphsOperImpl(expTree, set_tree_ownership)) 
{
    setTreeOperPimpl(pimpl);
}

GenerateMorphsOper::GenerateMorphsOper(bool set_tree_ownership) : 
pimpl(new GenerateMorphsOper::GenerateMorphsOperImpl(set_tree_ownership))
{
    setTreeOperPimpl(pimpl);
}

void GenerateMorphsOper::operator()() {
    (*pimpl)();
}

GenerateMorphsOper::GenerateMorphsOperImpl::GenerateMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree, bool set_tree_ownership) :
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, mSetTreeOwnershipForMorphs(set_tree_ownership)
{
    // no action
}

GenerateMorphsOper::GenerateMorphsOperImpl::GenerateMorphsOperImpl(bool set_tree_ownership) : 
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, mSetTreeOwnershipForMorphs(set_tree_ownership)
{
    // no action
}

void GenerateMorphsOper::GenerateMorphsOperImpl::CollectMorphs::MorphCollector(std::shared_ptr<MolpherMol> morph, void *functor) {
    CollectMorphs *collect =
            (CollectMorphs *) functor;
    (*collect)(morph);
}

GenerateMorphsOper::GenerateMorphsOperImpl::CollectMorphs::CollectMorphs(
    ConcurrentMolVector &morphs
    , std::shared_ptr<ExplorationTree> tree
    , bool set_ownership) 
:
mMorphs(morphs), 
mTree(tree),
mSetTreeOwnership(set_ownership)
{
    mCollectAttemptCount = 0;
}

void GenerateMorphsOper::GenerateMorphsOperImpl::CollectMorphs::operator()(std::shared_ptr<MolpherMol> morph) {
    ++mCollectAttemptCount; // atomic
    ConcurrentSmileSet::const_accessor dummy;
    if (mDuplicateChecker.insert(dummy, morph->getSMILES())) {
        if (mSetTreeOwnership) {
            morph->setOwner(mTree);
        }
        mMorphs.push_back(morph);
    } else {
        // ignore duplicate
    }
}

unsigned int GenerateMorphsOper::GenerateMorphsOperImpl::CollectMorphs::WithdrawCollectAttemptCount() {
    unsigned int ret = mCollectAttemptCount;
    mCollectAttemptCount = 0;
    return ret;
}

void GenerateMorphsOper::GenerateMorphsOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_pimpl = tree->pimpl;
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (tree_pimpl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree_pimpl->threadCnt);
        }

        ConcurrentMolVector leaves;
        tree_pimpl->fetchLeaves(tree, true, leaves);
        
        tree_pimpl->candidates.clear();
        CollectMorphs collectMorphs(tree_pimpl->candidates, tree, mSetTreeOwnershipForMorphs);
        for (auto leaf : leaves) {
            unsigned int morphAttempts = tree_pimpl->params.cntMorphs;
            if (leaf->getDistToTarget() < tree_pimpl->params.distToTargetDepthSwitch) {
                morphAttempts = tree_pimpl->params.cntMorphsInDepth;
            }

            tree_pimpl->candidates.reserve(tree_pimpl->candidates.size() + morphAttempts);

            std::vector<MolpherMol> decoys_dummy;
            std::vector<ChemOperSelector> oper_selectors;
            for (auto oper : tree_pimpl->chemOpers) {
                oper_selectors.push_back(static_cast<ChemOperSelector>(oper));
            }
            GenerateMorphs(
                    *leaf,
                    morphAttempts,
                    static_cast<FingerprintSelector>(tree_pimpl->fingerprint),
                    static_cast<SimCoeffSelector>(tree_pimpl->simCoeff),
                    oper_selectors,
                    tree_pimpl->target,
                    decoys_dummy,
                    tbbCtx,
                    &collectMorphs,
                    CollectMorphs::MorphCollector);
            MorphDerivationMap::accessor ac;

            if (tree_pimpl->morphDerivations.find(ac, leaf->getSMILES())) {
                ac->second += collectMorphs.WithdrawCollectAttemptCount();
            } else {
                tree_pimpl->morphDerivations.insert(ac, leaf->getSMILES());
                ac->second = collectMorphs.WithdrawCollectAttemptCount();
            }
        }
        tree_pimpl->candidates.shrink_to_fit();
        tree_pimpl->candidatesMask.clear();
        tree_pimpl->candidatesMask.resize(tree_pimpl->candidates.size(), true);
        SynchCout("Generated " + parseNumber(tree_pimpl->candidates.size()) + " morphs.");
    } else {
        throw std::runtime_error("Cannot generate morphs. No tree associated with this instance.");
    }
}

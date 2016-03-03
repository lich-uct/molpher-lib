
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "core/misc/inout.h"

#include "operations/ExtendTreeOper.hpp"
#include "core/API/ExplorationTreeImpl.h"
#include "ExtendTreeOperImpl.hpp"
#include "TreeOperationImpl.hpp"

ExtendTreeOper::ExtendTreeOper(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new ExtendTreeOper::ExtendTreeOperImpl(expTree))
{
    setTreeOperPimpl(pimpl);
}

ExtendTreeOper::ExtendTreeOper() :
pimpl(new ExtendTreeOper::ExtendTreeOperImpl())
{
    setTreeOperPimpl(pimpl);
}

void ExtendTreeOper::operator()() {
    (*pimpl)();
}


// pimpl

ExtendTreeOper::ExtendTreeOperImpl::ExtendTreeOperImpl(std::shared_ptr<ExplorationTree> expTree) :
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
{
    // no action
}

ExtendTreeOper::ExtendTreeOperImpl::ExtendTreeOperImpl() :
TreeOperation::TreeOperationImpl::TreeOperationImpl()
{
    // no action
}

ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::AcceptMorphs(
        ConcurrentMolVector &morphs, std::vector<bool> &survivors
        , std::shared_ptr<ExplorationTree> tree
        , ConcurrentSmileSet &modifiedParents
        ) :
mMorphs(morphs),
mSurvivors(survivors),
mTree(tree),
mTreePimpl(tree->pimpl),
mModifiedParents(modifiedParents),
mSurvivorCount(0) {
    assert(mMorphs.size() == mSurvivors.size());
}

ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::AcceptMorphs(
        AcceptMorphs &toSplit, tbb::split
        ) :
mTreePimpl(toSplit.mTreePimpl),
mMorphs(toSplit.mMorphs),
mSurvivors(toSplit.mSurvivors),
mModifiedParents(toSplit.mModifiedParents),
mSurvivorCount(0) {
}

void ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::operator()(
        const tbb::blocked_range<size_t> &r, tbb::pre_scan_tag) {
    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        if (mSurvivors[idx]) {
            ++mSurvivorCount;
        }
    }
}

void ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::operator()(
        const tbb::blocked_range<size_t> &r, tbb::final_scan_tag) {
    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        if (mSurvivors[idx]) {
            if (mSurvivorCount < mTreePimpl->params.cntCandidatesToKeepMax) {
                TreeMap::accessor ac;
                auto morph_smiles = mMorphs[idx]->getSMILES();

                if (mTreePimpl->treeMap.find(ac, morph_smiles)) {
                    SynchCout("Candidate morph: " + morph_smiles + " already present in the tree. Skipping...");
                    continue;
                }
            
                mTreePimpl->treeMap.insert(ac, morph_smiles);
                ac->second = mMorphs[idx];
                ac.release();

                if (mTreePimpl->treeMap.find(ac, mMorphs[idx]->getParentSMILES())) {
                    ac->second->addToDescendants(morph_smiles);
                    ac->second->addToHistoricDescendants(morph_smiles);
                    ConcurrentSmileSet::const_accessor dummy;
                    mModifiedParents.insert(dummy, ac->second->getSMILES());
                } else {
                    assert(false);
                }
            }

            ++mSurvivorCount;
        }
    }
}

void ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::reverse_join(AcceptMorphs &toJoin) {
    mSurvivorCount += toJoin.mSurvivorCount;
}

void ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs::assign(AcceptMorphs &toAssign) {
    mSurvivorCount = toAssign.mSurvivorCount;
}

ExtendTreeOper::ExtendTreeOperImpl::UpdateTree::UpdateTree(std::shared_ptr<ExplorationTree> tree) :
mTree(tree),
mTreePimpl(tree->pimpl)
{
}

void ExtendTreeOper::ExtendTreeOperImpl::UpdateTree::operator()(
        const ConcurrentSmileSet::range_type &modifiedParents) const {
    ConcurrentSmileSet::iterator itParent;
    for (itParent = modifiedParents.begin();
            itParent != modifiedParents.end(); itParent++) {

        // Determine what child is the closest to the target.
        double minDistance = DBL_MAX;
        TreeMap::accessor acParent;
        if (mTreePimpl->treeMap.find(acParent, itParent->first)) {

            std::set<std::string>::iterator itChild;
            auto descendants = acParent->second->getDescendants();
            for (itChild = descendants.begin();
                    itChild != descendants.end();
                    itChild++) {

                TreeMap::const_accessor acChild;
                if (mTreePimpl->treeMap.find(acChild, (*itChild))) {
                    double dist = acChild->second->getDistToTarget();
                    if (dist < minDistance) {
                        minDistance = dist;
                    }
                    if (!acChild->second->getTree()) {
                        acChild->second->setOwner(mTree);
                    } else {
                        assert(mTree == acChild->second->getTree());
                    }
                } else {
                    assert(false);
                }

            }

        } else {
            assert(false);
        }

        // Update the tree branch towards root.
        while (!acParent->second->getParentSMILES().empty()) {
            if (minDistance < acParent->second->getDistToTarget()) {
                acParent->second->setItersWithoutDistImprovement(0);
            }
            std::string smile = acParent->second->getParentSMILES();
            acParent.release();
            mTreePimpl->treeMap.find(acParent, smile);
            assert(!acParent.empty());
        }

    }
}

///**
// * Accept single morph with given index do not control anything.
// * @param idx Index or morph to accept.
// * @param morphs List of morphs.
// * @param ctx Context.
// * @param modifiedParents Parent to modify.
// */
//void ExtendTreeOper::ExtendTreeOperImpl::acceptMorph(
//        size_t idx,
//        ConcurrentMolVector &morphs,
//        std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_pimpl,
//        ConcurrentSmileSet &modifiedParents) {
//    auto morph_smiles = morphs[idx]->getSMILES();
//    TreeMap::accessor ac;
//    tree_pimpl->treeMap.insert(ac, );
//    ac->second = morphs[idx];
//    ac.release();
//
//    if (tree_pimpl->treeMap.find(ac, morphs[idx]->getSMILES())) {
//        ac->second->addToDescendants(morph_smiles);
//        ac->second->addToHistoricDescendants(morph_smiles);
//        ConcurrentSmileSet::const_accessor dummy;
//        modifiedParents.insert(dummy, ac->second->getSMILES());
//    } else {
//        assert(false);
//    }
//}

/**
 * Accept morphs from list. If there is no decoy the PathFinder::AcceptMorphs is 
 * used. Otherwise for each decoy the same number of best candidates is accepted.
 * @param morphs Candidates.
 * @param survivors Survive index.
 * @param ctx Context.
 * @param modifiedParents
 * @param decoySize Number of decoy used during exploration.
 */
void ExtendTreeOper::ExtendTreeOperImpl::acceptMorphs(ConcurrentMolVector &morphs,
        std::vector<bool> &survivors,
        std::shared_ptr<ExplorationTree> tree,
        ConcurrentSmileSet &modifiedParents,
        int decoySize) {

    // no decoy .. we can use old parallel approach        
    ExtendTreeOper::ExtendTreeOperImpl::AcceptMorphs acceptMorphs(morphs, survivors, tree, modifiedParents);
    // FIXME
    // Current TBB version does not support parallel_scan cancellation.
    // If it will be improved in the future, pass task_group_context
    // argument similarly as in parallel_for.
    tbb::parallel_scan(
            tbb::blocked_range<size_t>(0, morphs.size()),
            acceptMorphs, tbb::auto_partitioner());
    return;


    // Advance decoy functionality use if //if (decoySize != 0) {
    /*
    int maxDecoy = 0;
    for (size_t i = 0; i < morphs.size(); ++i) {
        if (survivors[i]) {
            maxDecoy = maxDecoy < morphs[i].nextDecoy ? morphs[i].nextDecoy : maxDecoy;
        }
    }
    // +1 is for target ..
    assert(maxDecoy < (decoySize + 1));
    // increase max decoy for target
    maxDecoy += 1;
    std::vector<int> limits(maxDecoy);
    // give the same to each
    std::stringstream ss; 
    ss << "number of active decoys: " << maxDecoy << std::endl;
    ss << "limits:      ";
    for (size_t i = 0; i < maxDecoy; ++i) {
        limits[i] = (ctx.params.cntCandidatesToKeepMax / maxDecoy) - 1;
        ss << limits[i] << " ";
    }
    ss << std::endl;
    
    bool print = true;
    // now scan the candidates list .. and store results
    for (size_t i = 0; i < morphs.size(); ++i) {
        if (survivors[i]) {
            if (print) {
                std::stringstream ss;
                if (decoySize < morphs[i].nextDecoy) {
                    ss << "next decoy:" << morphs[i].nextDecoy; 
                } else {
                    ss << "target";
                }
                ss <<  " distance:" << std::setprecision(5) << 
                        morphs[i].distToClosestDecoy;
                SynchCout(ss.str());
                print = false;
            }
            
            // determine the index to limits
            size_t limitIndex = morphs[i].nextDecoy;
            // do we have space left?
            if (limits[limitIndex] > 0) {
                --limits[limitIndex];
                acceptMorph(i, morphs, ctx, modifiedParents);
            } else {
                // we are out of space .. 
            }
        }
    }
    ss << "left unused: ";
    for (size_t i = 0; i < maxDecoy; ++i) {
        ss << limits[i] << " ";
    }
    SynchCout(ss.str());
     */
}

void ExtendTreeOper::ExtendTreeOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_pimpl = tree->pimpl;
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (tree_pimpl->threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(tree_pimpl->threadCnt);
        }

        ConcurrentSmileSet modifiedParents;
        acceptMorphs(
                tree_pimpl->candidates
                , tree_pimpl->candidatesMask
                , tree
                , modifiedParents
                , 0
                );

        UpdateTree updateTree(tree);
        tbb::parallel_for(ConcurrentSmileSet::range_type(modifiedParents),
                updateTree, tbb::auto_partitioner(), tbbCtx);

        tree_pimpl->candidates.clear();
        tree_pimpl->candidatesMask.clear();
        ConcurrentMolVector dummy;
        tree_pimpl->fetchLeaves(tree, true, dummy);
        tree_pimpl->generationCnt++;
    } else {
        throw std::runtime_error("Can't extend. No tree associated with this instance.");
    }
}



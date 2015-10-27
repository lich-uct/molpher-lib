
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>

#include "molpher_API/operations/ExtendTreeOper.hpp"

ExtendTreeOper::ExtendTreeOper(ExplorationTree& expTree) : TreeOperation(expTree) {
    // no action
}

ExtendTreeOper::AcceptMorphs::AcceptMorphs(
        ExplorationTree::MoleculeVector &morphs, std::vector<bool> &survivors,
        PathFinderContext &ctx, ExplorationTree::SmileSet &modifiedParents
        ) :
mMorphs(morphs),
mSurvivors(survivors),
mCtx(ctx),
mModifiedParents(modifiedParents),
mSurvivorCount(0) {
    assert(mMorphs.size() == mSurvivors.size());
}

ExtendTreeOper::AcceptMorphs::AcceptMorphs(
        AcceptMorphs &toSplit, tbb::split
        ) :
mCtx(toSplit.mCtx),
mMorphs(toSplit.mMorphs),
mSurvivors(toSplit.mSurvivors),
mModifiedParents(toSplit.mModifiedParents),
mSurvivorCount(0) {
}

void ExtendTreeOper::AcceptMorphs::operator()(
        const tbb::blocked_range<size_t> &r, tbb::pre_scan_tag) {
    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        if (mSurvivors[idx]) {
            ++mSurvivorCount;
        }
    }
}

void ExtendTreeOper::AcceptMorphs::operator()(
        const tbb::blocked_range<size_t> &r, tbb::final_scan_tag) {
    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        if (mSurvivors[idx]) {
            if (mSurvivorCount < mCtx.params.cntCandidatesToKeepMax) {
                PathFinderContext::CandidateMap::accessor ac;

                mCtx.candidates.insert(ac, mMorphs[idx].smile);
                ac->second = mMorphs[idx];
                ac.release();

                if (mCtx.candidates.find(ac, mMorphs[idx].parentSmile)) {
                    ac->second.descendants.insert(mMorphs[idx].smile);
                    ac->second.historicDescendants.insert(mMorphs[idx].smile);
                    ExplorationTree::SmileSet::const_accessor dummy;
                    mModifiedParents.insert(dummy, ac->second.smile);
                } else {
                    assert(false);
                }
            }

            ++mSurvivorCount;
        }
    }
}

void ExtendTreeOper::AcceptMorphs::reverse_join(AcceptMorphs &toJoin) {
    mSurvivorCount += toJoin.mSurvivorCount;
}

void ExtendTreeOper::AcceptMorphs::assign(AcceptMorphs &toAssign) {
    mSurvivorCount = toAssign.mSurvivorCount;
}

ExtendTreeOper::UpdateTree::UpdateTree(PathFinderContext &ctx) :
mCtx(ctx) {
}

void ExtendTreeOper::UpdateTree::operator()(
        const ExplorationTree::SmileSet::range_type &modifiedParents) const {
    ExplorationTree::SmileSet::iterator itParent;
    for (itParent = modifiedParents.begin();
            itParent != modifiedParents.end(); itParent++) {

        // Determine what child is the closest to the target.
        double minDistance = DBL_MAX;
        PathFinderContext::CandidateMap::accessor acParent;
        if (mCtx.candidates.find(acParent, itParent->first)) {

            std::set<std::string>::iterator itChild;
            for (itChild = acParent->second.descendants.begin();
                    itChild != acParent->second.descendants.end();
                    itChild++) {

                PathFinderContext::CandidateMap::const_accessor acChild;
                if (mCtx.candidates.find(acChild, (*itChild))) {
                    if (acChild->second.distToTarget < minDistance) {
                        minDistance = acChild->second.distToTarget;
                    }
                } else {
                    assert(false);
                }

            }

        } else {
            assert(false);
        }

        // Update the tree branch towards root.
        while (!acParent->second.parentSmile.empty()) {
            if (minDistance < acParent->second.distToTarget) {
                acParent->second.itersWithoutDistImprovement = 0;
            }
            std::string smile = acParent->second.parentSmile;
            acParent.release();
            mCtx.candidates.find(acParent, smile);
            assert(!acParent.empty());
        }

    }
}

/**
 * Accept single morph with given index do not control anything.
 * @param idx Index or morph to accept.
 * @param morphs List of morphs.
 * @param ctx Context.
 * @param modifiedParents Parent to modify.
 */
void ExtendTreeOper::acceptMorph(
        size_t idx,
        ExplorationTree::MoleculeVector &morphs,
        PathFinderContext &ctx,
        ExplorationTree::SmileSet &modifiedParents) {
    PathFinderContext::CandidateMap::accessor ac;
    ctx.candidates.insert(ac, morphs[idx].smile);
    ac->second = morphs[idx];
    ac.release();

    if (ctx.candidates.find(ac, morphs[idx].parentSmile)) {
        ac->second.descendants.insert(morphs[idx].smile);
        ac->second.historicDescendants.insert(morphs[idx].smile);
        ExplorationTree::SmileSet::const_accessor dummy;
        modifiedParents.insert(dummy, ac->second.smile);
    } else {
        assert(false);
    }
}

/**
 * Accept morphs from list. If there is no decoy the PathFinder::AcceptMorphs is 
 * used. Otherwise for each decoy the same number of best candidates is accepted.
 * @param morphs Candidates.
 * @param survivors Survive index.
 * @param ctx Context.
 * @param modifiedParents
 * @param decoySize Number of decoy used during exploration.
 */
void ExtendTreeOper::acceptMorphs(ExplorationTree::MoleculeVector &morphs,
        std::vector<bool> &survivors,
        PathFinderContext &ctx,
        ExplorationTree::SmileSet &modifiedParents,
        int decoySize) {

    // no decoy .. we can use old parallel approach        
    ExtendTreeOper::AcceptMorphs acceptMorphs(morphs, survivors, ctx, modifiedParents);
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

void ExtendTreeOper::operator()() {
    tbb::task_group_context tbbCtx;
    tbb::task_scheduler_init scheduler;
    if (threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(threadCnt);
    }

    PathFinderContext& context = fetchTreeContext();
    ExplorationTree::MoleculeVector& morphs = fetchGeneratedMorphs();
    ExplorationTree::BoolVector& survivors = fetchGeneratedMorphsMask();
    
    
    ExplorationTree::SmileSet modifiedParents;
    acceptMorphs(morphs, survivors, context, modifiedParents, context.decoys.size());

    UpdateTree updateTree(context);
    tbb::parallel_for(ExplorationTree::SmileSet::range_type(modifiedParents),
                updateTree, tbb::auto_partitioner(), tbbCtx);
    
    morphs.clear();
    survivors.clear();
}



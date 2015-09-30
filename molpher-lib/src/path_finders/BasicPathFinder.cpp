
#include <cassert>
#include <cmath>
#include <cfloat>
#include <string>
#include <sstream>

#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_exception.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"
#include "chem/morphing/Morphing.hpp"
#include "core/JobManager.h"
#include "path_finders/BasicPathFinder.hpp"

BasicPathFinder::BasicPathFinder(
        tbb::task_group_context *tbbCtx, JobManager *jobManager, int threadCnt
        ) :
mTbbCtx(tbbCtx),
mJobManager(jobManager),
mThreadCnt(threadCnt) {
}

BasicPathFinder::~BasicPathFinder() {
}

bool BasicPathFinder::Cancelled() {
    return mTbbCtx->is_group_execution_cancelled();
}

BasicPathFinder::FindLeaves::FindLeaves(MoleculeVector &leaves) :
mLeaves(leaves) {
}

void BasicPathFinder::FindLeaves::operator()(
        const PathFinderContext::CandidateMap::range_type &candidates) const {
    PathFinderContext::CandidateMap::iterator it;
    for (it = candidates.begin(); it != candidates.end(); it++) {
        if (!it->second.parentSmile.empty()) {
            it->second.itersWithoutDistImprovement++;
        }
        bool isLeaf = it->second.descendants.empty();
        if (isLeaf) {
            mLeaves.push_back(it->second);
        }
    }
}

BasicPathFinder::CollectMorphs::CollectMorphs(MoleculeVector &morphs) :
mMorphs(morphs) {
    mCollectAttemptCount = 0;
}

void BasicPathFinder::CollectMorphs::operator()(const MolpherMolecule &morph) {
    ++mCollectAttemptCount; // atomic
    SmileSet::const_accessor dummy;
    if (mDuplicateChecker.insert(dummy, morph.smile)) {
        mMorphs.push_back(morph);
    } else {
        // ignore duplicate
    }
}

unsigned int BasicPathFinder::CollectMorphs::WithdrawCollectAttemptCount() {
    unsigned int ret = mCollectAttemptCount;
    mCollectAttemptCount = 0;
    return ret;
}

void MorphCollector(MolpherMolecule *morph, void *functor) {
    BasicPathFinder::CollectMorphs *collect =
            (BasicPathFinder::CollectMorphs *) functor;
    (*collect)(*morph);
}

// return true if "a" is closes to target then "b"

bool BasicPathFinder::CompareMorphs::operator()(
        const MolpherMolecule &a, const MolpherMolecule &b) const {
    /* Morphs are rated according to their proximity to the connecting line
     between their closest decoy and the target (i.e. sum of both distances is
     minimal on the connecting line between decoy and target). When sums for
     both morphs are equal, it is possible (but not necessary) that both
     morphs lie on the same connecting line. In that case, morphs are
     rated only according to their proximity to the target. Such comparison
     should allow convergence to the target even in the late stages of the
     algorithm when majority of morphs lie on the connecting line between
     decoy closest to the target and the target itself. */

    double aSum = a.distToTarget + a.distToClosestDecoy;
    double bSum = b.distToTarget + b.distToClosestDecoy;

    bool approximatelyEqual = (
            fabs(aSum - bSum) <= (32 * DBL_EPSILON * fmax(fabs(aSum), fabs(bSum))));

    if (approximatelyEqual) {
        return a.distToTarget < b.distToTarget;
    } else {
        return aSum < bSum;
    }

    /**
     * Just sort based on distance to decoy .. or target. We prefer those
     * with target.
     */

    /* Experimental evaluation
    if (a.nextDecoy == -1 && b.nextDecoy == -1) {
        // both go for target, so take just their distance
        return a.distToTarget < b.distToTarget;
    } else if (a.nextDecoy == -1) {
        // a goes for target while b not -> a is greater
        return true;
    } else if (b.nextDecoy == -1) {
        // b goes for target while a not -> b is greater
        return false;
    } else {
        // decide based on nextDecoy
        if (a.nextDecoy == b.nextDecoy) {
            // go for same decoy use distance to closes decoy
            return a.distToClosestDecoy < b.distToClosestDecoy;
        } else {
            // greater is better
            return a.nextDecoy > b.nextDecoy;
        }
    }*/
}

BasicPathFinder::FilterMorphs::FilterMorphs(PathFinderContext &ctx,
        size_t globalMorphCount, MoleculeVector &morphs, std::vector<bool> &survivors
        ) :
mCtx(ctx),
mGlobalMorphCount(globalMorphCount),
mMorphs(morphs),
mSurvivors(survivors) {
    assert(mMorphs.size() == mSurvivors.size());
}

void BasicPathFinder::FilterMorphs::operator()(const tbb::blocked_range<size_t> &r) const {
    RDKit::RWMol* reqSubstructure = NULL;

    for (size_t idx = r.begin(); idx != r.end(); ++idx) {

        double acceptProbability = 1.0;
        bool isTarget = (mMorphs[idx].smile == mCtx.target.smile);
        if (idx >= mCtx.params.cntCandidatesToKeep && !isTarget) {
            acceptProbability =
                    0.25 - (idx - mCtx.params.cntCandidatesToKeep) /
                    ((mGlobalMorphCount - mCtx.params.cntCandidatesToKeep) * 4.0);
        }

        bool mightSurvive =
                SynchRand::GetRandomNumber(0, 99) < (int) (acceptProbability * 100);
        if (mightSurvive) {
            bool isDead = false;
            bool badWeight = false;
            bool badSascore = false;
            bool alreadyInTree = false;
            bool alreadyTriedByParent = false;
            bool tooManyProducedMorphs = false;
            bool badSubstructure = false;

            // Tests are ordered according to their cost.
            // Added test for SAScore

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);
            if (!isDead) {
                badWeight =
                        (mMorphs[idx].molecularWeight <
                        mCtx.params.minAcceptableMolecularWeight) ||
                        (mMorphs[idx].molecularWeight >
                        mCtx.params.maxAcceptableMolecularWeight);
            }

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);

            if (!isDead) {
                if (mCtx.params.useSyntetizedFeasibility) {
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
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);
            if (!isDead) {
                PathFinderContext::CandidateMap::const_accessor ac;
                if (mCtx.candidates.find(ac, mMorphs[idx].smile)) {
                    alreadyInTree = true;
                }
            }

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);
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
            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);
            if (!isDead) {
                PathFinderContext::MorphDerivationMap::const_accessor ac;
                if (mCtx.morphDerivations.find(ac, mMorphs[idx].smile)) {
                    tooManyProducedMorphs =
                            (ac->second > mCtx.params.cntMaxMorphs);
                }
            }

            isDead = (badWeight || badSascore || alreadyInTree ||
                    alreadyTriedByParent || tooManyProducedMorphs || badSubstructure);
            mSurvivors[idx] = !isDead;
        }
    }
    // release data
    if (reqSubstructure != NULL) {
        delete reqSubstructure;
    }
}

BasicPathFinder::AcceptMorphs::AcceptMorphs(
        MoleculeVector &morphs, std::vector<bool> &survivors,
        PathFinderContext &ctx, SmileSet &modifiedParents
        ) :
mMorphs(morphs),
mSurvivors(survivors),
mCtx(ctx),
mModifiedParents(modifiedParents),
mSurvivorCount(0) {
    assert(mMorphs.size() == mSurvivors.size());
}

BasicPathFinder::AcceptMorphs::AcceptMorphs(
        AcceptMorphs &toSplit, tbb::split
        ) :
mCtx(toSplit.mCtx),
mMorphs(toSplit.mMorphs),
mSurvivors(toSplit.mSurvivors),
mModifiedParents(toSplit.mModifiedParents),
mSurvivorCount(0) {
}

void BasicPathFinder::AcceptMorphs::operator()(
        const tbb::blocked_range<size_t> &r, tbb::pre_scan_tag) {
    for (size_t idx = r.begin(); idx != r.end(); ++idx) {
        if (mSurvivors[idx]) {
            ++mSurvivorCount;
        }
    }
}

void BasicPathFinder::AcceptMorphs::operator()(
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
                    SmileSet::const_accessor dummy;
                    mModifiedParents.insert(dummy, ac->second.smile);
                } else {
                    assert(false);
                }
            }

            ++mSurvivorCount;
        }
    }
}

void BasicPathFinder::AcceptMorphs::reverse_join(AcceptMorphs &toJoin) {
    mSurvivorCount += toJoin.mSurvivorCount;
}

void BasicPathFinder::AcceptMorphs::assign(AcceptMorphs &toAssign) {
    mSurvivorCount = toAssign.mSurvivorCount;
}

BasicPathFinder::UpdateTree::UpdateTree(PathFinderContext &ctx) :
mCtx(ctx) {
}

void BasicPathFinder::UpdateTree::operator()(
        const SmileSet::range_type &modifiedParents) const {
    BasicPathFinder::SmileSet::iterator itParent;
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

BasicPathFinder::PruneTree::PruneTree(PathFinderContext &ctx, SmileSet &deferred) :
mCtx(ctx),
mDeferred(deferred) {
}

void BasicPathFinder::PruneTree::operator()(
        const std::string &smile, tbb::parallel_do_feeder<std::string> &feeder) const {
    PathFinderContext::CandidateMap::accessor ac;
    mCtx.candidates.find(ac, smile);
    assert(!ac.empty());

    SmileSet::const_accessor dummy;
    bool deferred = mDeferred.find(dummy, smile);
    bool prune = (deferred ||
            (ac->second.itersWithoutDistImprovement > mCtx.params.itThreshold));
    if (prune) {

        bool tooManyDerivations = false;
        PathFinderContext::MorphDerivationMap::const_accessor acDerivations;
        if (mCtx.morphDerivations.find(acDerivations, smile)) {
            tooManyDerivations = (acDerivations->second > mCtx.params.cntMaxMorphs);
        }

        bool pruneThis = (deferred || tooManyDerivations);

        if (pruneThis) {
            PathFinderContext::CandidateMap::accessor acParent;
            mCtx.candidates.find(acParent, ac->second.parentSmile);
            assert(!acParent.empty());

            acParent->second.descendants.erase(smile);
            acParent.release();
            ac.release();

            EraseSubTree(smile);
        } else {
            std::set<std::string>::const_iterator it;
            for (it = ac->second.descendants.begin();
                    it != ac->second.descendants.end(); it++) {
                EraseSubTree(*it);
            }
            ac->second.descendants.clear();
            ac->second.itersWithoutDistImprovement = 0;
        }

    } else {
        std::set<std::string>::const_iterator it;
        for (it = ac->second.descendants.begin();
                it != ac->second.descendants.end(); it++) {
            feeder.add(*it);
        }
    }
}

void BasicPathFinder::PruneTree::EraseSubTree(const std::string &root) const {
    std::deque<std::string> toErase;
    toErase.push_back(root);

    while (!toErase.empty()) {
        std::string current = toErase.front();
        toErase.pop_front();

        PathFinderContext::CandidateMap::accessor ac;
        mCtx.candidates.find(ac, current);
        assert(!ac.empty());

        std::set<std::string>::const_iterator it;
        for (it = ac->second.descendants.begin();
                it != ac->second.descendants.end(); it++) {
            toErase.push_back(*it);
        }

        mCtx.prunedDuringThisIter.push_back(current);
        mCtx.candidates.erase(ac);
    }
}

BasicPathFinder::AccumulateTime::AccumulateTime(PathFinderContext &ctx) :
mCtx(ctx) {
    mTimestamp = std::clock();
}

unsigned int BasicPathFinder::AccumulateTime::GetElapsedSeconds(bool reset) {
    clock_t current = std::clock();
    unsigned int seconds =
            (unsigned int) ((current - mTimestamp) / CLOCKS_PER_SEC);
    if (reset) {
        mTimestamp = current;
    }
    return seconds;
}

void BasicPathFinder::AccumulateTime::ReportElapsedMiliseconds(
        const std::string &consumer, bool reset) {
    clock_t current = std::clock();
#if PATHFINDER_REPORTING == 1
    std::ostringstream stream;
    stream << mCtx.jobId << "/" << mCtx.iterIdx + 1 << ": " <<
            consumer << " consumed " << current - mTimestamp << " msec.";
    SynchCout(stream.str());
#endif
    if (reset) {
        mTimestamp = current;
    }
}

void BasicPathFinder::AccumulateTime::Reset() {
    mTimestamp = std::clock();
}

/**
 * Accept single morph with given index do not control anything.
 * @param idx Index or morph to accept.
 * @param morphs List of morphs.
 * @param ctx Context.
 * @param modifiedParents Parent to modify.
 */
void acceptMorph(
        size_t idx,
        BasicPathFinder::MoleculeVector &morphs,
        PathFinderContext &ctx,
        BasicPathFinder::SmileSet &modifiedParents) {
    PathFinderContext::CandidateMap::accessor ac;
    ctx.candidates.insert(ac, morphs[idx].smile);
    ac->second = morphs[idx];
    ac.release();

    if (ctx.candidates.find(ac, morphs[idx].parentSmile)) {
        ac->second.descendants.insert(morphs[idx].smile);
        ac->second.historicDescendants.insert(morphs[idx].smile);
        BasicPathFinder::SmileSet::const_accessor dummy;
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
void acceptMorphs(BasicPathFinder::MoleculeVector &morphs,
        std::vector<bool> &survivors,
        PathFinderContext &ctx,
        BasicPathFinder::SmileSet &modifiedParents,
        int decoySize) {

    // no decoy .. we can use old parallel approach        
    BasicPathFinder::AcceptMorphs acceptMorphs(morphs, survivors, ctx, modifiedParents);
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

void BasicPathFinder::operator()() {
    SynchCout(std::string("PathFinder thread started."));

    tbb::task_scheduler_init scheduler;
    if (mThreadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(mThreadCnt);
    }

    bool canContinueCurrentJob = false;
    bool pathFound = false;

    while (true) {

        if (!canContinueCurrentJob) {
            if (mJobManager->GetJob(mCtx)) {
                canContinueCurrentJob = true;
                pathFound = false;

                // Initialize the first iteration of a job.
                if (mCtx.candidates.empty()) {
                    assert(mCtx.iterIdx == 0);
                    PathFinderContext::CandidateMap::accessor ac;
                    mCtx.candidates.insert(ac, mCtx.source.smile);
                    ac->second = mCtx.source;
                }
            } else {
                break; // Thread termination.
            }
        }

        try {

            if (!Cancelled()) {
                mJobManager->GetFingerprintSelector(mCtx.fingerprintSelector);
                mJobManager->GetSimCoeffSelector(mCtx.simCoeffSelector);
                mJobManager->GetDimRedSelector(mCtx.dimRedSelector);
                mJobManager->GetChemOperSelectors(mCtx.chemOperSelectors);
                mJobManager->GetParams(mCtx.params);
                mJobManager->GetDecoys(mCtx.decoys);
                mCtx.prunedDuringThisIter.clear();
            }

            AccumulateTime molpherStopwatch(mCtx);
            AccumulateTime stageStopwatch(mCtx);

            MoleculeVector leaves;
            FindLeaves findLeaves(leaves);
            if (!Cancelled()) {
                tbb::parallel_for(
                        PathFinderContext::CandidateMap::range_type(mCtx.candidates),
                        findLeaves, tbb::auto_partitioner(), *mTbbCtx);
                stageStopwatch.ReportElapsedMiliseconds("FindLeaves", true);
            }

            /* TODO MPI
             MASTER
             prepare light snapshot (PathFinderContext::ContextToLightSnapshot)
             broadcast light snapshot
             convert leaves to std::vector
             scatter leaves over cluster
             convert node-specific part back to MoleculeVector

             SLAVE
             broadcast light snapshot
             convert light snapshot into context (PathFinderContext::SnapshotToContext)
             scatter leaves over cluster
             convert node-specific part back to MoleculeVector
             */

            MoleculeVector morphs;
            CollectMorphs collectMorphs(morphs);
            for (MoleculeVector::iterator it = leaves.begin(); it != leaves.end(); it++) {
                MolpherMolecule &candidate = (*it);
                unsigned int morphAttempts = mCtx.params.cntMorphs;
                if (candidate.distToTarget < mCtx.params.distToTargetDepthSwitch) {
                    morphAttempts = mCtx.params.cntMorphsInDepth;
                }

                if (!Cancelled()) {
                    morphs.reserve(morphs.size() + morphAttempts);

                    GenerateMorphs(
                            candidate,
                            morphAttempts,
                            mCtx.fingerprintSelector,
                            mCtx.simCoeffSelector,
                            mCtx.chemOperSelectors,
                            mCtx.target,
                            mCtx.decoys,
                            *mTbbCtx,
                            &collectMorphs,
                            MorphCollector);
                    PathFinderContext::MorphDerivationMap::accessor ac;

                    if (mCtx.morphDerivations.find(ac, candidate.smile)) {
                        ac->second += collectMorphs.WithdrawCollectAttemptCount();
                    } else {
                        mCtx.morphDerivations.insert(ac, candidate.smile);
                        ac->second = collectMorphs.WithdrawCollectAttemptCount();
                    }
                }

                if (Cancelled()) {
                    break;
                }
            }
            morphs.shrink_to_fit();

            if (!Cancelled()) {
                stageStopwatch.ReportElapsedMiliseconds("GenerateMorphs", true);
            }

            CompareMorphs compareMorphs;
            if (!Cancelled()) {
                /* FIXME
                 Current TBB version does not support parallel_sort cancellation.
                 If it will be improved in the future, pass task_group_context
                 argument similarly as in parallel_for. */
                tbb::parallel_sort(
                        morphs.begin(), morphs.end(), compareMorphs);
                stageStopwatch.ReportElapsedMiliseconds("SortMorphs", true);
            }

            /* TODO MPI
             MASTER
             gather snapshots from slaves
             update mCtx.morphDerivations according to gathered snapshots (consider using TBB for this)
             convert morphs to std::vector
             gather other morphs from slaves
               - each vector is pre-sorted and without duplicates
             integrate all morph vectors into final vector (consider using TBB for this)
               - check for cross-vector duplicates
               - merge sort

             SLAVE
             convert context to full snapshot (PathFinderContext::ContextToSnapshot)
             gather snapshot to master
             convert morphs to std::vector
             gather morphs back to master
             */

            /* TODO MPI
             MASTER
             prepare full snapshot (PathFinderContext::ContextToSnapshot)
             broadcast snapshot
             broadcast morph vector complete size
             scatter morph vector over cluster
             convert node-specific part back to MoleculeVector

             SLAVE
             broadcast snapshot
             convert snapshot into context (PathFinderContext::SnapshotToContext)
             broadcast morph vector complete size
             scatter morph vector over cluster
             convert node-specific part back to MoleculeVector
             */

            std::vector<bool> survivors;
            survivors.resize(morphs.size(), false);
            FilterMorphs filterMorphs(mCtx, morphs.size(), morphs, survivors);
            if (!Cancelled()) {
                if (mCtx.params.useSyntetizedFeasibility) {
                    SynchCout("\tUsing syntetize feasibility");
                }
                tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, morphs.size()),
                        filterMorphs, tbb::auto_partitioner(), *mTbbCtx);
                stageStopwatch.ReportElapsedMiliseconds("FilterMorphs", true);
            }

            /* TODO MPI
             MASTER
             gather survivors vector from slaves

             SLAVE
             gather survivors vector back to master
             */

            // Now we need to accept morphs ie. move the lucky one from 
            // morphs -> survivors
            SmileSet modifiedParents;
            acceptMorphs(morphs, survivors, mCtx, modifiedParents, mCtx.decoys.size());
            stageStopwatch.ReportElapsedMiliseconds("AcceptMorphs", true);

            UpdateTree updateTree(mCtx);
            if (!Cancelled()) {
                tbb::parallel_for(SmileSet::range_type(modifiedParents),
                        updateTree, tbb::auto_partitioner(), *mTbbCtx);
                stageStopwatch.ReportElapsedMiliseconds("UpdateTree", true);
            }

            if (!Cancelled()) {
                PathFinderContext::CandidateMap::const_accessor acTarget;
                mCtx.candidates.find(acTarget, mCtx.target.smile);
                pathFound = !acTarget.empty();
            }

            SmileSet deferredSmiles;
            SmileVector pruningQueue;
            PruneTree pruneTree(mCtx, deferredSmiles);
            if (!pathFound && !Cancelled()) {
                // Prepare deferred visual pruning.
                std::vector<MolpherMolecule> deferredMols;
                mJobManager->GetPruned(deferredMols);
                std::vector<MolpherMolecule>::iterator it;
                for (it = deferredMols.begin(); it != deferredMols.end(); it++) {
                    SmileSet::const_accessor dummy;
                    if (it->smile == mCtx.source.smile) {
                        continue;
                    }
                    deferredSmiles.insert(dummy, it->smile);
                }
                deferredMols.clear();

                pruningQueue.push_back(mCtx.source.smile);
                tbb::parallel_do(
                        pruningQueue.begin(), pruningQueue.end(), pruneTree, *mTbbCtx);
                stageStopwatch.ReportElapsedMiliseconds("PruneTree", true);
            }

            // find closest candidate
            PathFinderContext::CandidateMap::iterator itCandidates;
            double distance = DBL_MAX;
            for (itCandidates = mCtx.candidates.begin();
                    itCandidates != mCtx.candidates.end(); itCandidates++) {
                if (itCandidates->second.distToTarget < distance) {
                    distance = itCandidates->second.distToTarget;
                }
            }
            std::stringstream ss;
            ss << mCtx.jobId << "/" << mCtx.iterIdx << ": "
                    << "Closest candidate: " << morphs[0].distToTarget
                    << " -- " << morphs[0].smile << " -- ";
            SynchCout(ss.str());

            // find out if we should continue
            if (!Cancelled()) {
                mCtx.iterIdx += 1;
                mCtx.elapsedSeconds += molpherStopwatch.GetElapsedSeconds();

                if (canContinueCurrentJob) {
                    bool itersDepleted = (mCtx.params.cntIterations <= mCtx.iterIdx);
                    bool timeDepleted = (mCtx.params.timeMaxSeconds <= mCtx.elapsedSeconds);
                    canContinueCurrentJob = (!itersDepleted && !timeDepleted);
                }
            }
        } catch (tbb::tbb_exception &exc) {
            SynchCout(std::string(exc.what()));
            canContinueCurrentJob = false;
        }

        canContinueCurrentJob = mJobManager->CommitIteration(
                mCtx, canContinueCurrentJob, pathFound);

        SynchCout(std::string("Iteration commited."));

    }

    SynchCout(std::string("PathFinder thread terminated."));
}

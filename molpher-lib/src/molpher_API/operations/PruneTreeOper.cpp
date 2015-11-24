
#include "inout.h"

#include "molpher_API/operations/PruneTreeOper.hpp"
#include "molpher_API/callbacks/EraseSubtreeCallback.hpp"

PruneTreeOper::PruneTreeOper(ExplorationTree& expTree) : TreeOperation(expTree) {
    // no action
}

PruneTreeOper::PruneTree::PruneTree(PathFinderContext& ctx, TraverseCallback& callback) :
mCtx(ctx)
, mCallback(callback)
{
    // no action
}

void PruneTreeOper::PruneTree::operator()(const std::string& smile, tbb::parallel_do_feeder<std::string>& feeder) const {
    PathFinderContext::CandidateMap::accessor ac;
    mCtx.candidates.find(ac, smile);
    assert(!ac.empty());

    bool prune = (ac->second.itersWithoutDistImprovement > mCtx.params.itThreshold);
    if (prune) {

        bool tooManyDerivations = false;
        PathFinderContext::MorphDerivationMap::const_accessor acDerivations;
        if (mCtx.morphDerivations.find(acDerivations, smile)) {
            tooManyDerivations = (acDerivations->second > mCtx.params.cntMaxMorphs);
        }

        bool pruneThis = tooManyDerivations;

        if (pruneThis) {
            PathFinderContext::CandidateMap::accessor acParent;
            mCtx.candidates.find(acParent, ac->second.parentSmile);
            assert(!acParent.empty());

            acParent->second.descendants.erase(smile);
            acParent.release();
            ac.release();
            
            std::stringstream ss;
            ss << "Pruned: " << smile;
            SynchCout(ss.str());

            eraseSubtree(smile);
        } else {
            std::stringstream ss;
            ss << "Pruned (descendants only): " << smile;
            SynchCout(ss.str());
            
            std::set<std::string>::const_iterator it;
            for (it = ac->second.descendants.begin();
                    it != ac->second.descendants.end(); it++) {
                eraseSubtree(*it);
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

void PruneTreeOper::PruneTree::eraseSubtree(const std::string& smile) const {
    PathFinderContext::CandidateMap::accessor ac;
    mCtx.candidates.find(ac, smile);
    MolpherMol mol(ac->second);
    ac.release();
    mCallback.processMorph(mol);
}


void PruneTreeOper::operator()() {
    tbb::task_group_context tbbCtx;
    tbb::task_scheduler_init scheduler;
    if (threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(threadCnt);
    }
    ExplorationTree::SmileVector queue;
    queue.push_back(fetchTreeContext().source.smile);
    
    EraseSubtreeCallback callback(fetchTreeContext());
    PruneTree functor(fetchTreeContext(), callback);
    tbb::parallel_do(queue.begin(), queue.end(), functor, tbbCtx);
}

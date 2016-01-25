
#include "operations/TraverseOper.hpp"

TraverseOper::TraverseOper(ExplorationTree& expTree, TraverseCallback& callback) : TreeOperation(expTree), callback(callback), root(&(fetchTreeContext().source)) {
    // no action
}

TraverseOper::TraverseOper(TraverseCallback& callback) : TreeOperation(), callback(callback), root(nullptr) {
    // no action
}

TraverseOper::TraverseOper(ExplorationTree& expTree, TraverseCallback& callback, MolpherMolecule& root) : TreeOperation(expTree), callback(callback), root(&root)  {
    // no action
}

TraverseOper::TraverseOper(ExplorationTree& expTree, TraverseCallback& callback, MolpherMol& root) : TraverseOper(expTree, callback, root.fetchMolpherMolecule())  {
    // no action
}

TraverseOper::TraversalFunctor::TraversalFunctor(PathFinderContext& ctx, TraverseCallback& callback) : 
mCtx(ctx)
, mCallback(callback)
{
    // no action
}

void TraverseOper::TraversalFunctor::operator()(const std::string& smile, tbb::parallel_do_feeder<std::string>& feeder) const {
        PathFinderContext::CandidateMap::accessor ac;
        mCtx.candidates.find(ac, smile);
        assert(!ac.empty());
        
        makeCallback(ac->second);
        
        std::set<std::string>::const_iterator it;
        for (it = ac->second.descendants.begin();
                it != ac->second.descendants.end(); it++) {
            feeder.add(*it);
        }
}

void TraverseOper::TraversalFunctor::makeCallback(MolpherMolecule& morph) const {
    MolpherMol mol(morph);
    mCallback.processMorph(mol);
}

void TraverseOper::operator()() {
    if (this->tree) {
        tbb::task_group_context tbbCtx;
        tbb::task_scheduler_init scheduler;
        if (threadCnt > 0) {
            scheduler.terminate();
            scheduler.initialize(threadCnt);
        }
        ExplorationTree::SmileVector queue;
        if (!root) {
            root = &(fetchTreeContext().source);
        }
        queue.push_back(root->smile);

        TraversalFunctor functor(fetchTreeContext(), callback);
        tbb::parallel_do(queue.begin(), queue.end(), functor, tbbCtx);
    } else {
        throw std::runtime_error("Cannot traverse the tree. None associated with this instance.");
    }
}



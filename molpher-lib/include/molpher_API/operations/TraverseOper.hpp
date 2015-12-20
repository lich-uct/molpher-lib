/* 
 * File:   TraverseOper.hpp
 * Author: sichom
 *
 * Created on October 29, 2015, 9:46 AM
 */

#ifndef TRAVERSEOPER_HPP
#define	TRAVERSEOPER_HPP

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_do.h>

#include "TreeOperation.hpp"
#include "../callbacks/TraverseCallback.hpp"

class TraverseOper : public TreeOperation {
    
    private:
        class TraversalFunctor {
        public:
            TraversalFunctor(PathFinderContext &ctx, TraverseCallback& callback);
            void operator()(const std::string &smile,
                    tbb::parallel_do_feeder<std::string> &feeder) const;

        protected:
            void makeCallback(MolpherMolecule& root) const;

        private:
            PathFinderContext &mCtx;
            TraverseCallback &mCallback;
        };
        
        MolpherMolecule& root;
        TraverseCallback& callback;
        
        TraverseOper(ExplorationTree& expTree, TraverseCallback& callback, MolpherMolecule& root);
    
    public:
        TraverseOper(ExplorationTree& expTree, TraverseCallback& callback);
        TraverseOper(TraverseCallback& callback);
        TraverseOper(ExplorationTree& expTree, TraverseCallback& callback, MolpherMol& root);
        virtual void operator()();
};


#endif	/* TRAVERSEOPER_HPP */


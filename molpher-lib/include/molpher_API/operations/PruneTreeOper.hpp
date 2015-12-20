/* 
 * File:   PruneTreeOper.hpp
 * Author: sichom
 *
 * Created on October 29, 2015, 2:05 PM
 */

#ifndef PRUNETREEOPER_HPP
#define	PRUNETREEOPER_HPP

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_do.h>

#include "TreeOperation.hpp"
#include "../callbacks/TraverseCallback.hpp"

class PruneTreeOper : public TreeOperation {
    
    private:
        class PruneTree {
        public:
            PruneTree(PathFinderContext &ctx, TraverseCallback& callback);
            void operator()(const std::string &smile,
                    tbb::parallel_do_feeder<std::string> &feeder) const;

        protected:
            void eraseSubtree(const std::string& smile) const;

        private:
            PathFinderContext &mCtx;
            TraverseCallback &mCallback;
        };
    
    public:
        PruneTreeOper(ExplorationTree& expTree);
        PruneTreeOper();
        virtual void operator()();
};

#endif	/* PRUNETREEOPER_HPP */


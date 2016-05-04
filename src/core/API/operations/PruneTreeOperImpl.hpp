/* 
 * File:   PruneTreeOperImpl.hpp
 * Author: sichom
 *
 * Created on March 4, 2016, 4:51 PM
 */

#ifndef PRUNETREEOPERIMPL_HPP
#define	PRUNETREEOPERIMPL_HPP

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_do.h>

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class PruneTreeOper::PruneTreeOperImpl : public TreeOperation::TreeOperationImpl {
    
    private:
        class PruneTree {
        public:
            PruneTree(std::shared_ptr<ExplorationTree> expTree);
            void operator()(const std::string &smile,
                    tbb::parallel_do_feeder<std::string> &feeder) const;

        protected:
            void eraseSubtree(const std::string& smile, bool descendents_only) const;

        private:
            std::shared_ptr<ExplorationTree> mTree;
        };
    
    public:
        PruneTreeOperImpl(std::shared_ptr<ExplorationTree> expTree);
        PruneTreeOperImpl();
        virtual void operator()();
};

#endif	/* PRUNETREEOPERIMPL_HPP */


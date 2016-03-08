/* 
 * File:   TraverseOperImpl.hpp
 * Author: sichom
 *
 * Created on March 8, 2016, 11:35 AM
 */

#ifndef TRAVERSEOPERIMPL_HPP
#define	TRAVERSEOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class TraverseOper::TraverseOperImpl : public TreeOperation::TreeOperationImpl {
    
    private:
        class TraversalFunctor {
        public:
            TraversalFunctor(std::shared_ptr<ExplorationTree> expTree, std::shared_ptr<TraverseCallback> callback);
            void operator()(const std::string &smile,
                    tbb::parallel_do_feeder<std::string> &feeder) const;

        protected:
            void makeCallback(std::shared_ptr<MolpherMol> root) const;

        private:
            std::shared_ptr<ExplorationTree> mTree;
            std::shared_ptr<TraverseCallback> mCallback;
        };
        
        std::shared_ptr<MolpherMol> root;
        std::shared_ptr<TraverseCallback> callback;
    
    public:
        TraverseOperImpl(std::shared_ptr<ExplorationTree> expTree, std::shared_ptr<TraverseCallback> callback);
        TraverseOperImpl(std::shared_ptr<TraverseCallback> callback);
        TraverseOperImpl(std::shared_ptr<ExplorationTree> expTree, std::shared_ptr<TraverseCallback> callback, const std::string& rootSMILES);
        virtual void operator()();
};

#endif	/* TRAVERSEOPERIMPL_HPP */


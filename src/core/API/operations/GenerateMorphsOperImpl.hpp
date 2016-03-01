/* 
 * File:   GenerateMorphsOperImpl.hpp
 * Author: sichom
 *
 * Created on March 1, 2016, 9:06 AM
 */

#ifndef GENERATEMORPHSOPERIMPL_HPP
#define	GENERATEMORPHSOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class GenerateMorphsOper::GenerateMorphsOperImpl : public TreeOperation::TreeOperationImpl {

private:
    
    class CollectMorphs
    {
    public:
        CollectMorphs(ConcurrentMolVector &morphs, std::shared_ptr<ExplorationTree> tree);
        void operator()(std::shared_ptr<MolpherMol> morph);
        unsigned int WithdrawCollectAttemptCount();
        static void MorphCollector(std::shared_ptr<MolpherMol> morph, void *functor);

    private:
        ConcurrentSmileSet mDuplicateChecker;
        ConcurrentMolVector &mMorphs;
        std::shared_ptr<ExplorationTree> mTree;
        tbb::atomic<unsigned int> mCollectAttemptCount;
    };
    
public:
    GenerateMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree);
    GenerateMorphsOperImpl();
    virtual void operator()();

};

#endif	/* GENERATEMORPHSOPERIMPL_HPP */


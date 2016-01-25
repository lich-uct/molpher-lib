/* 
 * File:   PutativeExtend.hpp
 * Author: sichom
 *
 * Created on October 19, 2015, 8:39 AM
 */

#ifndef PUTATIVEEXTEND_HPP
#define	PUTATIVEEXTEND_HPP

#include "TreeOperation.hpp"

class GenerateMorphsOper : public TreeOperation {

private:
    
    class CollectMorphs
    {
    public:
        CollectMorphs(ExplorationTree::MoleculeVector &morphs);
        void operator()(const MolpherMolecule &morph);
        unsigned int WithdrawCollectAttemptCount();
        static void MorphCollector(MolpherMolecule *morph, void *functor);

    private:
        ExplorationTree::SmileSet mDuplicateChecker;
        ExplorationTree::MoleculeVector &mMorphs;
        tbb::atomic<unsigned int> mCollectAttemptCount;
    };
    
public:
    GenerateMorphsOper(ExplorationTree& expTree);
    GenerateMorphsOper();
    virtual void operator()();

};

#endif	/* PUTATIVEEXTEND_HPP */


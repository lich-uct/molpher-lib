/* 
 * File:   ExtendTreeOperImpl.hpp
 * Author: sichom
 *
 * Created on March 2, 2016, 2:57 PM
 */

#ifndef EXTENDTREEOPERIMPL_HPP
#define	EXTENDTREEOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class ExtendTreeOper::ExtendTreeOperImpl : public TreeOperation::TreeOperationImpl {
    
private:
    
    class AcceptMorphs
    {
    public:
        AcceptMorphs(ConcurrentMolVector &morphs, std::vector<bool> &survivors,
            std::shared_ptr<ExplorationTree> tree, ConcurrentSmileSet &modifiedParents);
        AcceptMorphs(AcceptMorphs &toSplit, tbb::split);
        void operator()(const tbb::blocked_range<size_t> &r, tbb::pre_scan_tag);
        void operator()(const tbb::blocked_range<size_t> &r, tbb::final_scan_tag);
        void reverse_join(AcceptMorphs &toJoin);
        void assign(AcceptMorphs &toAssign);

    private:
        ConcurrentMolVector &mMorphs;
        std::vector<bool> &mSurvivors;
        std::shared_ptr<ExplorationTree> mTree;
        std::shared_ptr<ExplorationTree::ExplorationTreeImpl> mTreePimpl;
        ConcurrentSmileSet &mModifiedParents;
        unsigned int mSurvivorCount;
    };
    
    class UpdateTree
    {
    public:
        UpdateTree(std::shared_ptr<ExplorationTree> tree);
        void operator()(const ConcurrentSmileSet::range_type &modifiedParents) const;

    private:
        std::shared_ptr<ExplorationTree> mTree;
        std::shared_ptr<ExplorationTree::ExplorationTreeImpl> mTreePimpl;
    };
    
//    static void acceptMorph(
//        size_t idx,
//        ConcurrentMolVector &morphs,
//        std::shared_ptr<ExplorationTree::ExplorationTreeImpl>,
//        ConcurrentSmileSet &modifiedParents);
    
    static void acceptMorphs(ConcurrentMolVector &morphs,
        std::vector<bool> &survivors,
        std::shared_ptr<ExplorationTree> tree,
        ConcurrentSmileSet &modifiedParents,
        int decoySize);
    
public:
    ExtendTreeOperImpl(std::shared_ptr<ExplorationTree> expTree);
    ExtendTreeOperImpl();
    virtual void operator()();
};

#endif	/* EXTENDTREEOPERIMPL_HPP */


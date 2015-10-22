/* 
 * File:   ExtendTreeOper.hpp
 * Author: sichom
 *
 * Created on October 22, 2015, 5:02 PM
 */

#ifndef EXTENDTREEOPER_HPP
#define	EXTENDTREEOPER_HPP

#include <tbb/parallel_scan.h>

#include "TreeOperation.hpp"

class ExtendTreeOper : public TreeOperation {
    
private:
    
    class AcceptMorphs
    {
    public:
        AcceptMorphs(ExplorationTree::MoleculeVector &morphs, std::vector<bool> &survivors,
            PathFinderContext &ctx, ExplorationTree::SmileSet &modifiedParents);
        AcceptMorphs(AcceptMorphs &toSplit, tbb::split);
        void operator()(const tbb::blocked_range<size_t> &r, tbb::pre_scan_tag);
        void operator()(const tbb::blocked_range<size_t> &r, tbb::final_scan_tag);
        void reverse_join(AcceptMorphs &toJoin);
        void assign(AcceptMorphs &toAssign);

    private:
        ExplorationTree::MoleculeVector &mMorphs;
        std::vector<bool> &mSurvivors;
        PathFinderContext &mCtx;
        ExplorationTree::SmileSet &mModifiedParents;
        unsigned int mSurvivorCount;
    };
    
    class UpdateTree
    {
    public:
        UpdateTree(PathFinderContext &ctx);
        void operator()(const ExplorationTree::SmileSet::range_type &modifiedParents) const;

    private:
        PathFinderContext &mCtx;
    };
    
    static void acceptMorph(
        size_t idx,
        ExplorationTree::MoleculeVector &morphs,
        PathFinderContext &ctx,
        ExplorationTree::SmileSet &modifiedParents);
    
    static void acceptMorphs(ExplorationTree::MoleculeVector &morphs,
        std::vector<bool> &survivors,
        PathFinderContext &ctx,
        ExplorationTree::SmileSet &modifiedParents,
        int decoySize);
    
public:
    ExtendTreeOper(ExplorationTree& expTree);
    virtual void operator()();
};

#endif	/* EXTENDTREEOPER_HPP */


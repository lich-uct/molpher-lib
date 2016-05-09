/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EXTENDTREEOPERIMPL_HPP
#define	EXTENDTREEOPERIMPL_HPP

#include <tbb/parallel_scan.h>

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


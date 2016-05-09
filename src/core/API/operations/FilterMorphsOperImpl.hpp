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

#ifndef FILTERMORPHSOPERIMPL_HPP
#define	FILTERMORPHSOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"

class FilterMorphsOper::FilterMorphsOperImpl : public TreeOperation::TreeOperationImpl {
private:
    int filters;
    bool verbose;

    class FilterMorphs {
    public:
        FilterMorphs(std::shared_ptr<ExplorationTree::ExplorationTreeImpl> tree_pimpl, size_t globalMorphCount,
                ConcurrentMolVector &morphs, std::vector<bool> &survivors, int filters, bool verbose);
        void operator()(const tbb::blocked_range<size_t> &r) const;

    private:
        std::shared_ptr<ExplorationTree::ExplorationTreeImpl> mTreePimpl;
        size_t mGlobalMorphCount;
        ConcurrentMolVector &mMorphs;
        std::vector<bool> &mSurvivors;
        int mFilters;
        bool mVerboseOutput;
    };

public:

    enum MorphFilters : int {
        PROBABILITY = 1 << 0,
        WEIGHT = 1 << 1,
        SYNTHESIS = 1 << 2,
        MAX_DERIVATIONS = 1 << 3,
        DUPLICATES = 1 << 4,
        HISTORIC_DESCENDENTS = 1 << 5,
        ALL = PROBABILITY | WEIGHT | SYNTHESIS | MAX_DERIVATIONS | DUPLICATES | HISTORIC_DESCENDENTS
    };

    FilterMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree, bool verbose = false);
    FilterMorphsOperImpl(bool verbose = false);
    FilterMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree, int filters, bool verbose = false);
    FilterMorphsOperImpl(int filters, bool verbose = false);
    virtual void operator()();
};

#endif	/* FILTERMORPHSOPERIMPL_HPP */


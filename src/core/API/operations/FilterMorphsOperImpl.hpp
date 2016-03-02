/* 
 * File:   FilterMorphsOperImpl.hpp
 * Author: sichom
 *
 * Created on March 2, 2016, 12:46 PM
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


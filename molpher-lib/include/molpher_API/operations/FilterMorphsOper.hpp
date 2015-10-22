/* 
 * File:   FilterMorphsOper.hpp
 * Author: sichom
 *
 * Created on October 22, 2015, 10:18 AM
 */

#ifndef FILTERMORPHSOPER_HPP
#define	FILTERMORPHSOPER_HPP

#include "TreeOperation.hpp"

class FilterMoprhsOper : public TreeOperation {
private:
    int filters;

    class FilterMorphs {
    public:
        FilterMorphs(PathFinderContext &ctx, size_t globalMorphCount,
                ExplorationTree::MoleculeVector &morphs, std::vector<bool> &survivors, int filters);
        void operator()(const tbb::blocked_range<size_t> &r) const;

    private:
        PathFinderContext &mCtx;
        size_t mGlobalMorphCount;
        ExplorationTree::MoleculeVector &mMorphs;
        std::vector<bool> &mSurvivors;
        int mFilters;
    };

public:

    enum MorphFilters : int {
        PROBABILITY = 1 << 0,
        WEIGHT = 1 << 1,
        SYNTHESIS = 1 << 2,
        COUNT = 1 << 3,
        ALL = PROBABILITY | WEIGHT | SYNTHESIS | COUNT
    };

    FilterMoprhsOper(ExplorationTree& expTree);
    FilterMoprhsOper(ExplorationTree& expTree, int filters);
    virtual void operator()();
};


#endif	/* FILTERMORPHSOPER_HPP */


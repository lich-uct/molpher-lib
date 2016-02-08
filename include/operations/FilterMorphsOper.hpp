///* 
// * File:   FilterMorphsOper.hpp
// * Author: sichom
// *
// * Created on October 22, 2015, 10:18 AM
// */
//
//#ifndef FILTERMORPHSOPER_HPP
//#define	FILTERMORPHSOPER_HPP
//
//#include "TreeOperation.hpp"
//
//class FilterMorphsOper : public TreeOperation {
//private:
//    int filters;
//    bool verbose;
//
//    class FilterMorphs {
//    public:
//        FilterMorphs(PathFinderContext &ctx, size_t globalMorphCount,
//                ExplorationTree::MoleculeVector &morphs, std::vector<bool> &survivors, int filters, bool verbose);
//        void operator()(const tbb::blocked_range<size_t> &r) const;
//
//    private:
//        PathFinderContext &mCtx;
//        size_t mGlobalMorphCount;
//        ExplorationTree::MoleculeVector &mMorphs;
//        std::vector<bool> &mSurvivors;
//        int mFilters;
//        bool mVerboseOutput;
//    };
//
//public:
//
//    enum MorphFilters : int {
//        PROBABILITY = 1 << 0,
//        WEIGHT = 1 << 1,
//        SYNTHESIS = 1 << 2,
//        MAX_DERIVATIONS = 1 << 3,
//        DUPLICATES = 1 << 4,
//        HISTORIC_DESCENDENTS = 1 << 5,
//        ALL = PROBABILITY | WEIGHT | SYNTHESIS | MAX_DERIVATIONS | DUPLICATES | HISTORIC_DESCENDENTS
//    };
//
//    FilterMorphsOper(ExplorationTree& expTree, bool verbose = false);
//    FilterMorphsOper(bool verbose = false);
//    FilterMorphsOper(ExplorationTree& expTree, int filters, bool verbose = false);
//    FilterMorphsOper(int filters, bool verbose = false);
//    virtual void operator()();
//};
//
//
//#endif	/* FILTERMORPHSOPER_HPP */
//

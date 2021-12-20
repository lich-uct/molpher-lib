/*
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

#pragma once

#include <string>
#include <vector>
#include <memory>

#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>
#include <morphing/operators/MorphingOperator.hpp>

#include "core/misc/global_types.h"

#include "core/data_structs/MolpherParam.h"
#include "data_structs/ExplorationTree.hpp"
#include "operations/TreeOperation.hpp"
#include "data_structs/MolpherMol.hpp"
#include "operations/FindLeavesOper.hpp"
#include "operations/GenerateMorphsOper.hpp"
#include "operations/SortMorphsOper.hpp"
#include "operations/FilterMorphsOper.hpp"
#include "operations/ExtendTreeOper.hpp"
#include "operations/PruneTreeOper.hpp"
#include "operations/TraverseOper.hpp"
#include "operations/CleanMorphsOper.hpp"

class ExplorationTree::ExplorationTreeImpl
{
    
    friend class TreeOperation::TreeOperationImpl;
    friend class FindLeavesOper::FindLeavesOperImpl;
    friend class GenerateMorphsOper::GenerateMorphsOperImpl;
    friend class SortMorphsOper::SortMorphsOperImpl;
    friend class FilterMorphsOper::FilterMorphsOperImpl;
    friend class ExtendTreeOper::ExtendTreeOperImpl;
    friend class PruneTreeOper::PruneTreeOperImpl;
    friend class TraverseOper::TraverseOperImpl;
    friend class CleanMorphsOper::CleanMorphsOperImpl;

    private:
        
        /**
         * Number of generations.
         */
        unsigned generationCnt;
        
        /**
         * Number of generations.
         */
        unsigned threadCnt;

        /**
         * Fingerprint used in algorithm.
         * @see FingerprintSelector
         */
        int fingerprint;

        /**
         * Similarity coef used in algorithm.
         * @see SimCoeffSelector
         */
        int simCoeff;

        /**
         * Vector with used morphing operators. Determines
         * possible morphing operation applied during generation
         * of new morphs.
         * @see MorphingOperator
         */
        std::vector<std::shared_ptr<MorphingOperator>> chemOpers;

        /**
         * Parameters for morphing algorithm.
         */
        MolpherParam params;

        /**
         * Source molecule.
         */
        std::shared_ptr<MolpherMol> source;

        /**
         * Target molecule.
         */
        std::shared_ptr<MolpherMol> target;

    //    typedef tbb::concurrent_vector<std::string> PrunedVector;

        /**
         * Candidate morphs.
         */
        ConcurrentMolVector candidates;

        /**
         * Candidates mask.
         */
		ConcurrentMaskVector candidatesMask;

        /**
         * Molecules in the tree.
         */
        TreeMap treeMap;

        MorphDerivationMap morphDerivations;

//        PrunedVector pruned;
        
        void erase(const std::string& canonSMILES);
        
    public:
        ExplorationTreeImpl();
        
        std::shared_ptr<ExplorationData> asData() const;
        void updateData(const ExplorationData& data, std::shared_ptr<ExplorationTree> tree);

        const std::vector<std::shared_ptr<MorphingOperator>>& getMorphingOperators();
        void setMorphingOperators(const std::vector<std::shared_ptr<MorphingOperator>>& operators);
        void addMorphingOperator(std::shared_ptr<MorphingOperator> operator_);
        
        void runOperation(TreeOperation& operation, std::shared_ptr<ExplorationTree> tree);
        MolVector fetchLeaves(std::shared_ptr<ExplorationTree> tree, bool increase_dist_improve_counter);
        void fetchLeaves(std::shared_ptr<ExplorationTree> tree, bool increase_dist_improve_counter, ConcurrentMolVector& ret);
        std::shared_ptr<MolpherMol> fetchMol(const std::string& canonSMILES);
        bool hasMol(const std::string& canonSMILES);
        bool hasMol(std::shared_ptr<MolpherMol> mol);
        bool isPathFound();
        void deleteSubtree(const std::string& canonSMILES, bool descendents_only);
        void generateMorphs(std::shared_ptr<ExplorationTree>);
        void generateMorphs(std::shared_ptr<ExplorationTree>, const std::vector<std::shared_ptr<MorphCollector> >& collectors);
        void sortMorphs(std::shared_ptr<ExplorationTree>);
        void filterMorphs(std::shared_ptr<ExplorationTree> tree, bool verbose_output);
        void filterMorphs(FilterMorphsOper::MorphFilters filters, std::shared_ptr<ExplorationTree> tree, bool verbose_output);
        void extend(std::shared_ptr<ExplorationTree> tree);
        void prune(std::shared_ptr<ExplorationTree> tree);
        void traverse(std::shared_ptr<ExplorationTree> tree, const std::string& rootSMILES, TraverseCallback& callback);
        void traverse(std::shared_ptr<ExplorationTree> tree, TraverseCallback& callback);
        void save(const std::string& filename);
        
        int getThreadCount();
        MolVector getCandidateMorphs();
        std::vector<bool> getCandidateMorphsMask(); // TODO add a bitset version
        unsigned getGenerationCount();
        
        void setThreadCount(int threadCnt);
        void setCandidateMorphsMask(const std::vector<bool>&); // TODO add a bitset version
};

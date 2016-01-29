/*
 Copyright (c) 2012 Petr Koupy

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

#include "global_types.h"
//#include "core/misc/selectors/fingerprint_selectors.h"
//#include "core/misc/selectors/simcoeff_selectors.h"
//#include "dimred_selectors.h"
//#include "core/misc/selectors/chemoper_selectors.h"

#include "data_structs/MolpherParam.h"
#include "core/data_structs/ExplorationData.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "operations/TreeOperation.hpp"
#include "data_structs/MolpherMol.hpp"

//class TreeOperation::TreeOperationImpl; // forward declaration to resolve circular dependency

class ExplorationTree::ExplorationTreeImpl
{
    
    friend class TreeOperation::TreeOperationImpl;

    private:
        
        /**
         * Number of generations.
         */
        unsigned generationCnt;

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
         * Vector or used morphing operator. Determine
         * possible morphing operation applied during generating
         * new morphs.
         * @see ChemOperSelector
         */
        std::vector<int> chemOpers;

        /**
         * Parameters for morphing algorithm.
         */
        MolpherParam params;

        /**
         * Source molecule.
         */
        std::shared_ptr<MolpherMol::MolpherMolImpl> source;

        /**
         * Target molecule.
         */
        std::shared_ptr<MolpherMol::MolpherMolImpl> target;

    //    typedef tbb::concurrent_vector<std::string> PrunedVector;

        /**
         * Candidate morphs.
         */
        ConcurrentMolVector candidates;

        /**
         * Candidates mask.
         */
        std::vector<bool> candidatesMask;

        /**
         * Molecules in the tree.
         */
        TreeMap treeMap;

        MorphDerivationMap morphDerivations;

//        PrunedVector pruned;
        
    public:
        
        static std::shared_ptr<ExplorationTree::ExplorationTreeImpl> createFromData(const ExplorationData &data);
        
        ExplorationTreeImpl(const std::string& sourceMolAsSMILES);
        ExplorationTreeImpl(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES);
        ExplorationTreeImpl(const ExplorationData &data);
        
        std::shared_ptr<ExplorationData> asData();
        updateFromData(const ExplorationData &data);
        
        void runOperation(std::shared_ptr<TreeOperation::TreeOperationImpl> operation);
        
        void fetchLeaves(std::vector<std::shared_ptr<MolpherMol::MolpherMolImpl> >& leaves, bool increase_dist_improve_counter = false);
        std::shared_ptr<const std::vector<std::shared_ptr<MolpherMol::MolpherMolImpl> > > fetchLeaves();
        std::shared_ptr<MolpherMol::MolpherMolImpl> fetchMol(const std::string& canonSMILES);
        bool hasMol(const std::string& canonSMILES);
        bool isPathFound();
        void deleteSubtree(const std::string& canonSMILES);
        void generateMorphs();
        void sortMorphs();
        void filterMorphs();
        void filterMorphs(int filters);
        void extend();
        void prune();
        
        std::shared_ptr<const std::vector<std::shared_ptr<MolpherMol::MolpherMolImpl> > > getCandidateMorphs();
        const std::vector<bool>& getCandidateMorphsMask(); // TODO add a bitset version
        
        void setCandidateMorphsMask(const std::vector<bool>&); // TODO add a bitset version
};

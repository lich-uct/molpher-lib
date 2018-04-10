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

#ifndef EXPLORATIONDATAIMPL_HPP
#define	EXPLORATIONDATAIMPL_HPP

#pragma once

#include <boost/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>
#include <core/misc/inout.h>

#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"
#include "selectors/chemoper_selectors.h"

#include "core/misc/global_types.h"

#include "MolpherParam.h"
#include "data_structs/ExplorationData.hpp"
#include "MolpherMolData.hpp"

struct ExplorationData::ExplorationDataImpl
{
//    typedef std::vector<std::string> PrunedVector;
    
    /**
     * Number of generations.
     */
    unsigned generationCnt;
    
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
     * Vector or used morphing operator. Determine
     * possible morphing operation applied during generating
     * new morphs.
     * @see ChemOperSelector
     */
    std::set<int> chemOpers;

    /**
     * Parameters for morphing algorithm.
     */
    MolpherParam params;

    /**
     * Source molecule.
     */
    MolpherMolData source;
    
    /**
     * Target molecule.
     */
    MolpherMolData target;
    
    /**
     * Candidate morphs.
     */
    CandidatesVectorData candidates;
    
    /**
     * Candidates mask.
     */
    CandidatesMaskVectorData candidatesMask;

    /**
     * Molecules in the tree.
     */
    TreeMapData treeMap;
    
    MorphDerivationMapData morphDerivations;
    
//    PrunedVector pruned;

    friend class boost::serialization::access;
    
    bool isValid() const {
        if (chemOpers.empty()) {
            Cerr("No chemical operators specified.");
            return false;
        }
        if (!params.isValid()) {
            Cerr("Parameters instance is not valid.");
            return false;
        }
        if (!source.isValid()) {
            Cerr("Source is not a valid instance.");
            return false;
        }
//        if (!target.isValid()) {
//            Cerr("Target is not a valid instance.");
//            return false;
//        }
        if (treeMap.empty()) {
            Cerr("There are no molecules in the tree.");
            return false;
        }
        if (source.SMILES.compare(target.SMILES) == 0) {
            Cerr("Source and target are the same.");
            return false;
        }
        return true;
    }

    ExplorationDataImpl() :
            generationCnt(0)
            , threadCnt(0)
        {
            fingerprint = DEFAULT_FP;
            simCoeff = DEFAULT_SC;

            // add all operators by default
            chemOpers.insert(OP_ADD_ATOM);
            chemOpers.insert(OP_REMOVE_ATOM);
            chemOpers.insert(OP_ADD_BOND);
            chemOpers.insert(OP_REMOVE_BOND);
            chemOpers.insert(OP_MUTATE_ATOM);
            chemOpers.insert(OP_INTERLAY_ATOM);
            chemOpers.insert(OP_BOND_REROUTE);
            chemOpers.insert(OP_BOND_CONTRACTION);
        }

    template <class Archive>
    void save(Archive & ar, const unsigned int version) const {

        ar  & BOOST_SERIALIZATION_NVP(generationCnt)
            & BOOST_SERIALIZATION_NVP(threadCnt) 
            & BOOST_SERIALIZATION_NVP(fingerprint) 
            & BOOST_SERIALIZATION_NVP(simCoeff) 
            & BOOST_SERIALIZATION_NVP(chemOpers) 
            & BOOST_SERIALIZATION_NVP(params) 
            & BOOST_SERIALIZATION_NVP(source) 
            & BOOST_SERIALIZATION_NVP(target) 
            & BOOST_SERIALIZATION_NVP(treeMap) 
            & BOOST_SERIALIZATION_NVP(candidates)
            & BOOST_SERIALIZATION_NVP(candidatesMask)
            & BOOST_SERIALIZATION_NVP(morphDerivations);
    //        & BOOST_SERIALIZATION_NVP(pruned);
    }

    template <class Archive>
    void load(Archive & ar, const unsigned int version) {

        ar  & BOOST_SERIALIZATION_NVP(generationCnt)
            & BOOST_SERIALIZATION_NVP(threadCnt) 
            & BOOST_SERIALIZATION_NVP(fingerprint) 
            & BOOST_SERIALIZATION_NVP(simCoeff) 
            & BOOST_SERIALIZATION_NVP(chemOpers) 
            & BOOST_SERIALIZATION_NVP(params) 
            & BOOST_SERIALIZATION_NVP(source) 
            & BOOST_SERIALIZATION_NVP(target) 
            & BOOST_SERIALIZATION_NVP(treeMap) 
            & BOOST_SERIALIZATION_NVP(candidates)
            & BOOST_SERIALIZATION_NVP(candidatesMask)
            & BOOST_SERIALIZATION_NVP(morphDerivations);
    //        & BOOST_SERIALIZATION_NVP(pruned);

    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

// turn off versioning
BOOST_CLASS_IMPLEMENTATION(ExplorationData, object_serializable)
// turn off tracking
BOOST_CLASS_TRACKING(ExplorationData, track_never)

#endif	/* EXPLORATIONDATAIMPL_HPP */


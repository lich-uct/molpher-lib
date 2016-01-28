/* 
 * File:   ExplorationData.hpp
 * Author: sichom
 *
 * Created on January 27, 2016, 10:11 AM
 */

#ifndef EXPLORATIONDATA_HPP
#define	EXPLORATIONDATA_HPP

#pragma once

#include <string>
#include <vector>
#include <map>

#include <boost/version.hpp>
//#include <boost/cstdint.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include "core/misc/selectors/fingerprint_selectors.h"
#include "core/misc/selectors/simcoeff_selectors.h"
#include "core/misc/selectors/chemoper_selectors.h"

#include "MolpherParam.h"
#include "MolpherMolData.hpp"

struct ExplorationData
{
    typedef std::vector<MolpherMolData> CandidatesVector;
    typedef std::vector<bool> CandidatesMaskVector;
    typedef std::map<std::string, MolpherMolData> TreeMap;
    typedef std::map<std::string, unsigned> MorphDerivationMap;
//    typedef std::vector<std::string> PrunedVector;
    
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
    MolpherMolData source;
    
    /**
     * Target molecule.
     */
    MolpherMolData target;
    
    /**
     * Candidate morphs.
     */
    CandidatesVector candidates;
    
    /**
     * Candidates mask.
     */
    CandidatesMaskVector candidatesMask;

    /**
     * Molecules in the tree.
     */
    TreeMap treeMap;
    
    MorphDerivationMap morphDerivations;
    
//    PrunedVector pruned;
    
    ExplorationData();
    bool isValid();

    friend class boost::serialization::access;
    
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;
    template<class Archive>
    void load(Archive & ar, const unsigned int version);
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

bool ExplorationData::isValid() {
    return (!chemOpers.empty()) &&
            params.isValid() &&
            source.isValid() &&
            target.isValid() &&
            (source.SMILES != target.SMILES);
}

ExplorationData::ExplorationData() :
        generationCnt(0)
    {
        fingerprint = DEFAULT_FP;
        simCoeff = DEFAULT_SC;
        
        // add all operators by default
        chemOpers.push_back(OP_ADD_ATOM);
        chemOpers.push_back(OP_REMOVE_ATOM);
        chemOpers.push_back(OP_ADD_BOND);
        chemOpers.push_back(OP_REMOVE_BOND);
        chemOpers.push_back(OP_MUTATE_ATOM);
        chemOpers.push_back(OP_INTERLAY_ATOM);
        chemOpers.push_back(OP_BOND_REROUTE);
        chemOpers.push_back(OP_BOND_CONTRACTION);
    }

// turn off versioning
BOOST_CLASS_IMPLEMENTATION(ExplorationData, object_serializable)
// turn off tracking
BOOST_CLASS_TRACKING(ExplorationData, track_never)

template <class Archive>
void ExplorationData::save(Archive & ar, const unsigned int version) const {
    
    ar  & BOOST_SERIALIZATION_NVP(generationCnt) 
        & BOOST_SERIALIZATION_NVP(fingerprint) 
        & BOOST_SERIALIZATION_NVP(simCoeff) 
        & BOOST_SERIALIZATION_NVP(chemOpers) 
        & BOOST_SERIALIZATION_NVP(params) 
        & BOOST_SERIALIZATION_NVP(source) 
        & BOOST_SERIALIZATION_NVP(target) 
        & BOOST_SERIALIZATION_NVP(treeMap) 
        & BOOST_SERIALIZATION_NVP(morphDerivations) 
        & BOOST_SERIALIZATION_NVP(pruned);
}

template <class Archive>
void ExplorationData::load(Archive & ar, const unsigned int version) const {
    
    ar  & BOOST_SERIALIZATION_NVP(generationCnt) 
        & BOOST_SERIALIZATION_NVP(fingerprint) 
        & BOOST_SERIALIZATION_NVP(simCoeff) 
        & BOOST_SERIALIZATION_NVP(chemOpers) 
        & BOOST_SERIALIZATION_NVP(params) 
        & BOOST_SERIALIZATION_NVP(source) 
        & BOOST_SERIALIZATION_NVP(target) 
        & BOOST_SERIALIZATION_NVP(treeMap) 
        & BOOST_SERIALIZATION_NVP(morphDerivations) 
        & BOOST_SERIALIZATION_NVP(pruned);
    
}

#endif	/* EXPLORATIONDATA_HPP */


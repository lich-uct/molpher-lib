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
#include <map>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include "core/misc/selectors/fingerprint_selectors.h"
#include "core/misc/selectors/simcoeff_selectors.h"
//#include "dimred_selectors.h"
#include "core/misc/selectors/chemoper_selectors.h"

#include "MolpherParam.h"
#include "MolpherMolecule.h"

struct IterationSnapshot
{
    IterationSnapshot() :
        jobId(0),
        iterIdx(0),
        elapsedSeconds(0)
    {
        fingerprintSelector = DEFAULT_FP;
        simCoeffSelector = DEFAULT_SC;
//        dimRedSelector = DEFAULT_DR;
        // in default add every ChemOperSelector
        chemOperSelectors.push_back(OP_ADD_ATOM);
        chemOperSelectors.push_back(OP_REMOVE_ATOM);
        chemOperSelectors.push_back(OP_ADD_BOND);
        chemOperSelectors.push_back(OP_REMOVE_BOND);
        chemOperSelectors.push_back(OP_MUTATE_ATOM);
        chemOperSelectors.push_back(OP_INTERLAY_ATOM);
        chemOperSelectors.push_back(OP_BOND_REROUTE);
        chemOperSelectors.push_back(OP_BOND_CONTRACTION);
    }

    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        // usage of BOOST_SERIALIZATION_NVP macro enable us to use xml serialisation
        ar & BOOST_SERIALIZATION_NVP(jobId) &
            BOOST_SERIALIZATION_NVP(iterIdx) & 
            BOOST_SERIALIZATION_NVP(elapsedSeconds) &
            BOOST_SERIALIZATION_NVP(fingerprintSelector) & 
            BOOST_SERIALIZATION_NVP(simCoeffSelector) &
//            BOOST_SERIALIZATION_NVP(dimRedSelector) & 
            BOOST_SERIALIZATION_NVP(chemOperSelectors) &
            BOOST_SERIALIZATION_NVP(params) & 
            BOOST_SERIALIZATION_NVP(source) &
            BOOST_SERIALIZATION_NVP(target) & 
            BOOST_SERIALIZATION_NVP(decoys) &
            BOOST_SERIALIZATION_NVP(candidates) & 
            BOOST_SERIALIZATION_NVP(morphDerivations) &
            BOOST_SERIALIZATION_NVP(prunedDuringThisIter);
    }
    
    bool IsValid()
    {
        bool decoysValid = true;
        std::vector<MolpherMolecule>::iterator it;
        for (it = decoys.begin(); it != decoys.end(); ++it) {
            if (!it->IsValid()) {
                decoysValid = false;
                break;
            }
        }

        return (!chemOperSelectors.empty()) &&
            params.IsValid() &&
            source.IsValid() &&
            target.IsValid() &&
            (source.smile != target.smile) &&
            decoysValid;
    }

    typedef std::map<std::string, MolpherMolecule> CandidateMap;
    typedef std::map<std::string, boost::uint32_t> MorphDerivationMap;
    typedef std::vector<std::string> PrunedMoleculeVector;
    
    /**
     * Job id.
     */
    boost::uint32_t jobId;
    
    /**
     * Iteration id.
     */
    boost::uint32_t iterIdx;
    
    /**
     * Total time spend by working on this instance of snapshot.
     */
    boost::uint32_t elapsedSeconds;

    /**
     * Fingerprint used in algorithm.
     * @see FingerprintSelector
     */
    boost::int32_t fingerprintSelector;
    
    /**
     * Similarity coef used in algorithm.
     * @see SimCoeffSelector
     */
    boost::int32_t simCoeffSelector;
    
    /**
     * Id of used location computing algorithm.
     * @see DimRedSelector
     */
    boost::int32_t dimRedSelector;
    
    /**
     * Vector or used morphing operator. Determine
     * possible morphing operation applied during generating
     * new morphs.
     * @see ChemOperSelector
     */
    std::vector<boost::int32_t> chemOperSelectors;

    /**
     * Parameters for morphing algorithm.
     */
    MolpherParam params;

    /**
     * Source molecule.
     */
    MolpherMolecule source;
    
    /**
     * Target molecule.
     */
    MolpherMolecule target;
       
    /**
     * Decoys used during exploration of chemical space.
     */
    std::vector<MolpherMolecule> decoys;

    /**
     * Candidate molecules ie. molecule storage.
     */
    CandidateMap candidates;
    
    MorphDerivationMap morphDerivations;
    
    PrunedMoleculeVector prunedDuringThisIter;

};

// add information about version to archive
//BOOST_CLASS_IMPLEMENTATION(IterationSnapshot, object_class_info)
// turn off versioning
BOOST_CLASS_IMPLEMENTATION(IterationSnapshot, object_serializable)
// turn off tracking
BOOST_CLASS_TRACKING(IterationSnapshot, track_never)
// specify version
//BOOST_CLASS_VERSION(IterationSnapshot, 1)

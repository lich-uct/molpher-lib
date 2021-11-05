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

#include <boost/version.hpp>
#include <boost/cstdint.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/nvp.hpp>

struct MolpherParam
{
    MolpherParam() :
        cntCandidatesToKeep(50),
        cntCandidatesToKeepMax(100),
        cntMorphs(80),
        cntMorphsInDepth(150),
        distToTargetDepthSwitch(0.15),
        cntMaxMorphs(1500),
        itThreshold(5),
        cntIterations(1000),
        timeMaxSeconds(21600000),
        minAcceptableMolecularWeight(0.0),
        maxAcceptableMolecularWeight(100000.0),
        useSyntetizedFeasibility(true),
        useSubstructureRestriction(false),
        decoyRange(0.2),
		sascoreMax(6.0)
    {
    }

    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        // usage of BOOST_SERIALIZATION_NVP enable us to use xml serialisation
        ar & BOOST_SERIALIZATION_NVP(cntCandidatesToKeep);
        ar & BOOST_SERIALIZATION_NVP(cntCandidatesToKeepMax);
        ar & BOOST_SERIALIZATION_NVP(cntMorphs);
        ar & BOOST_SERIALIZATION_NVP(cntMorphsInDepth);
        ar & BOOST_SERIALIZATION_NVP(distToTargetDepthSwitch);
        ar & BOOST_SERIALIZATION_NVP(cntMaxMorphs);
        ar & BOOST_SERIALIZATION_NVP(itThreshold);
        ar & BOOST_SERIALIZATION_NVP(cntIterations);
        ar & BOOST_SERIALIZATION_NVP(timeMaxSeconds);
        ar & BOOST_SERIALIZATION_NVP(minAcceptableMolecularWeight);
        ar & BOOST_SERIALIZATION_NVP(maxAcceptableMolecularWeight);
        ar & BOOST_SERIALIZATION_NVP(useSyntetizedFeasibility);
        
        //if (version > 0) {
            ar & BOOST_SERIALIZATION_NVP(useSubstructureRestriction);
            ar & BOOST_SERIALIZATION_NVP(decoyRange);
		ar & BOOST_SERIALIZATION_NVP(sascoreMax);
        //}
    }

    bool isValid() const
    {
        return (cntCandidatesToKeepMax > 0) &&
            (cntCandidatesToKeepMax >= cntCandidatesToKeep) &&
            (cntMorphs > 0) &&
            (cntMorphsInDepth > 0) &&
            (distToTargetDepthSwitch >= 0.0) &&
            (distToTargetDepthSwitch <= 1.0) &&
            (cntMaxMorphs > 0) &&
            (itThreshold > 0) &&
            (cntIterations > 0) &&
            (timeMaxSeconds > 0) &&
            (minAcceptableMolecularWeight >= 0.0) &&
            (minAcceptableMolecularWeight <= maxAcceptableMolecularWeight) &&
            (maxAcceptableMolecularWeight > 0.0) &&
            (decoyRange >= 0) && (decoyRange <= 1.0) &&
			(sascoreMax > 0);
    }

    // Molpher BIBE2011 parameters (do not rename)
    boost::uint32_t cntCandidatesToKeep;
	boost::uint32_t cntCandidatesToKeepMax;
	boost::uint32_t cntMorphs;
	boost::uint32_t cntMorphsInDepth;
	double distToTargetDepthSwitch;
    boost::uint32_t cntMaxMorphs;
	boost::uint32_t itThreshold;
    boost::uint32_t cntIterations;
	boost::uint32_t timeMaxSeconds;
    
    // newly added parameters to improve Molpher behaviour
    double minAcceptableMolecularWeight;
    double maxAcceptableMolecularWeight;
    
    /**
     * Filter molecules based on syntetized feasibility?
     */
    bool useSyntetizedFeasibility;
    
    /**
     * Use substructure filter? Added in version 1.
     */
    bool useSubstructureRestriction;
    
    /**
     * Decoy range. Ie. how close must candidate be to the
     * decoy to get over it. Must be from interval (1,0).
     */
    double decoyRange;

	/**
	 * Maximum allowed value of SAScore (synthetic accessibility score)
	 */
	double sascoreMax;
};

// add information about version to archive
//BOOST_CLASS_IMPLEMENTATION(MolpherParam, object_class_info)
// turn off versioning
BOOST_CLASS_IMPLEMENTATION(MolpherParam, object_serializable)
// turn off tracking        
BOOST_CLASS_TRACKING(MolpherParam, track_never)
//  set version
//BOOST_CLASS_VERSION(MolpherParam, 1)
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
#include <set>
#include <cfloat>

#include <boost/cstdint.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include <iostream>

struct MolpherMolecule
{
    MolpherMolecule() :
        parentChemOper(0),
        distToTarget(DBL_MAX),
        distToClosestDecoy(0),
        molecularWeight(0.0),
        sascore(0.0),
        itersWithoutDistImprovement(0),
        posX(0),
        posY(0)
    {
    }

    MolpherMolecule(std::string &smile) :
        smile(smile),
        parentChemOper(0),
        distToTarget(DBL_MAX),
        distToClosestDecoy(0),
        molecularWeight(0.0),
        sascore(0.0),
        itersWithoutDistImprovement(0),
        posX(0),
        posY(0)
    {
    }

    MolpherMolecule(std::string &smile, std::string &formula) :
        smile(smile),
        formula(formula),
        parentChemOper(0),
        distToTarget(DBL_MAX),
        distToClosestDecoy(0),
        molecularWeight(0.0),
        sascore(0.0),
        itersWithoutDistImprovement(0),
        posX(0),
        posY(0)
    {
    }

    MolpherMolecule(std::string &smile, std::string &formula,
        std::string &parentSmile, boost::int32_t parentChemOper,
        double distToTarget, double distToClosestDecoy, double molecularWeight, double sascore
        ) :
        smile(smile),
        formula(formula),
        parentChemOper(parentChemOper),
        parentSmile(parentSmile),
        distToTarget(distToTarget),
        distToClosestDecoy(distToClosestDecoy),
        molecularWeight(molecularWeight),
        sascore(sascore),
        itersWithoutDistImprovement(0),
        posX(0),
        posY(0)
    {
    }
    
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        // usage of BOOST_SERIALIZATION_NVP enable us to use xml serialisation
        ar & BOOST_SERIALIZATION_NVP(smile) & 
            BOOST_SERIALIZATION_NVP(formula) & 
            BOOST_SERIALIZATION_NVP(parentChemOper) & 
            BOOST_SERIALIZATION_NVP(parentSmile) & 
            BOOST_SERIALIZATION_NVP(descendants) &
            BOOST_SERIALIZATION_NVP(historicDescendants) & 
            BOOST_SERIALIZATION_NVP(distToTarget) & 
            BOOST_SERIALIZATION_NVP(distToClosestDecoy) &
            BOOST_SERIALIZATION_NVP(molecularWeight) & 
            BOOST_SERIALIZATION_NVP(itersWithoutDistImprovement) & 
            BOOST_SERIALIZATION_NVP(posX) & 
            BOOST_SERIALIZATION_NVP(posY);
    }

    bool IsValid()
    {
        return (!smile.empty());
    }

    std::string smile;
    std::string formula;
    boost::int32_t parentChemOper;
    std::string parentSmile;
    std::set<std::string> descendants;
    std::set<std::string> historicDescendants;

    double distToTarget;
    
    /**
     * Hold distance to the nextDecoy. Value -1 say that 
     * there is decoy that we should visit.
     */
    double distToClosestDecoy;
    double molecularWeight;
    
    /**
     * Keep sascore for filter, computed in generate morphs, 
     * constructors updated. The value is not serialised.
     */
    double sascore;
    
    boost::uint32_t itersWithoutDistImprovement;

    double posX;
    double posY;
    
    /**
     * Store index of next decoy that should be visited by 
     * this molecule. If greater then number of decoys 
     * then molecule can go straight for the target.
     * 
     * The value is not serialised.
     */
    // int nextDecoy;
};

// turn off versioning
BOOST_CLASS_IMPLEMENTATION(MolpherMolecule, object_serializable)
// turn off tracking
BOOST_CLASS_TRACKING(MolpherMolecule, track_never)

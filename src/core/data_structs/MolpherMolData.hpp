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

#ifndef MOLPHERMOLDATA_HPP
#define	MOLPHERMOLDATA_HPP

#include <string>
#include <set>
#include <cfloat>

#include <boost/serialization/set.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include <iostream>

struct MolpherMolData
{
    
    std::string SMILES;
    std::string formula;
    double molecularWeight;
    unsigned parentOper;
    std::string parentSmile;
    
    /**
     * Direct descendants of this molecule.
     */
    std::set<std::string> descendants;
    
    /**
     * All morphs created from this molecule.
     */
    std::set<std::string> historicDescendants;
    
    /**
     *  Syntetic accessibility score.
     */
    double sascore;
    
    /**
     *  Distance to target.
     */
    double distToTarget; 
    
    /**
     * Number of new generations created from this molecule without 
     * an improvement in distance.
     */
    unsigned gensWithoutDistImprovement;
    
    friend class boost::serialization::access;

    MolpherMolData() :
            parentOper(0),
            distToTarget(DBL_MAX),
            molecularWeight(0.0),
            sascore(0.0),
            gensWithoutDistImprovement(0)
    {
        // no action
    }
    
    MolpherMolData(const std::string& smiles) :
    SMILES(smiles),
    parentOper(0),
    distToTarget(DBL_MAX),
    molecularWeight(0.0),
    sascore(0.0),
    gensWithoutDistImprovement(0)
    {
        // no action
    }

    bool isValid() const {
        return (!SMILES.empty());
    }

    template <class Archive>
    void save(Archive & ar, const unsigned int version) const {
        auto& descendants_saved = descendants;
        auto& historicDescendants_saved = descendants;
        
        ar  & BOOST_SERIALIZATION_NVP(SMILES)
            & BOOST_SERIALIZATION_NVP(formula)
            & BOOST_SERIALIZATION_NVP(sascore)
            & BOOST_SERIALIZATION_NVP(parentOper) 
            & BOOST_SERIALIZATION_NVP(parentSmile)
            & BOOST_SERIALIZATION_NVP(descendants_saved) 
            & BOOST_SERIALIZATION_NVP(historicDescendants_saved) 
            & BOOST_SERIALIZATION_NVP(distToTarget) 
            & BOOST_SERIALIZATION_NVP(molecularWeight) 
            & BOOST_SERIALIZATION_NVP(gensWithoutDistImprovement);
    }

    template <class Archive>
    void load(Archive & ar, const unsigned int version) {

        std::set<std::string> descendants_saved;
        std::set<std::string> historicDescendants_saved;

        ar  & BOOST_SERIALIZATION_NVP(SMILES) 
            & BOOST_SERIALIZATION_NVP(formula)
            & BOOST_SERIALIZATION_NVP(sascore)
            & BOOST_SERIALIZATION_NVP(parentOper) 
            & BOOST_SERIALIZATION_NVP(parentSmile)
            & BOOST_SERIALIZATION_NVP(descendants_saved) 
            & BOOST_SERIALIZATION_NVP(historicDescendants_saved) 
            & BOOST_SERIALIZATION_NVP(distToTarget)
            & BOOST_SERIALIZATION_NVP(molecularWeight) 
            & BOOST_SERIALIZATION_NVP(gensWithoutDistImprovement);

        descendants.swap(descendants_saved);
        historicDescendants.swap(historicDescendants_saved);

    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

// turn off versioning
BOOST_CLASS_IMPLEMENTATION(MolpherMolData, object_serializable)
// turn off tracking
BOOST_CLASS_TRACKING(MolpherMolData, track_never)

#endif	/* MOLPHERMOLDATA_HPP */


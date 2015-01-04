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
#include <ctime>

#include <boost/cstdint.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/time_serialize.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>

#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "dimred_selectors.h"
#include "chemoper_selectors.h"

#include "MolpherMolecule.h"

struct NeighborhoodTask
{
    NeighborhoodTask() :
        attemptCount(10),
        maxDepth(1),
        maxDistance(0.1)
    {
        fingerprintSelector = DEFAULT_FP;
        simCoeffSelector = DEFAULT_SC;
        dimRedSelector = DEFAULT_DR;

        // Should help frontend to match request against backend response.
        taskTimestamp = boost::posix_time::microsec_clock::universal_time();
    }

    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & taskTimestamp & fingerprintSelector & simCoeffSelector &
            dimRedSelector & chemOperSelectors & origin & context &
            attemptCount & maxDepth & maxDistance;
    }

    bool IsValid()
    {
        bool contextValid = true;
        std::vector<MolpherMolecule>::iterator it;
        for (it = context.begin(); it != context.end(); ++it) {
            if (!it->IsValid()) {
                contextValid = false;
                break;
            }
        }

        bool withOrigin = (!chemOperSelectors.empty()) &&
            origin.IsValid() &&
            contextValid &&
            (attemptCount > 0) &&
            (maxDepth > 0) &&
            (maxDistance >= 0.00) &&
            (maxDistance <= 1.00);

        bool withoutOrigin = origin.smile.empty() && contextValid;

        return withOrigin || withoutOrigin;
    }

    boost::posix_time::ptime taskTimestamp;

    boost::int32_t fingerprintSelector;
    boost::int32_t simCoeffSelector;
    boost::int32_t dimRedSelector;
    std::vector<boost::int32_t> chemOperSelectors;

    MolpherMolecule origin;
    std::vector<MolpherMolecule> context;
    boost::uint32_t attemptCount;
    boost::uint32_t maxDepth;
    double maxDistance;
};

BOOST_CLASS_IMPLEMENTATION(NeighborhoodTask, object_serializable) // turn off versioning
BOOST_CLASS_TRACKING(NeighborhoodTask, track_never) // turn off tracking

struct NeighborhoodTaskResult
{
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & taskTimestamp & origin & reducedNeighborhood & reducedContext;
    }

    boost::posix_time::ptime taskTimestamp;
    MolpherMolecule origin;
    std::vector<MolpherMolecule> reducedNeighborhood;
    std::vector<MolpherMolecule> reducedContext;
};

BOOST_CLASS_IMPLEMENTATION(NeighborhoodTaskResult, object_serializable) // turn off versioning
BOOST_CLASS_TRACKING(NeighborhoodTaskResult, track_never) // turn off tracking

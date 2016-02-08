/*
 Copyright (c) 2012 Peter Szepe

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


#include <vector>
#include <string>

#include <GraphMol/GraphMol.h>

#include <tbb/atomic.h>
#include <tbb/blocked_range.h>

#include "selectors/chemoper_selectors.h"
#include "core/chem/fingerprintStrategy/FingerprintStrategy.h"
#include "core/chem/simCoefStrategy/SimCoefStrategy.h"
#include "MorphingData.h"
#include "core/chem/morphingStrategy/MorphingStrategy.h"
#include "core/chem/SimCoefCalculator.hpp"

class CalculateMorphs
{
public:
    CalculateMorphs(
        MorphingData &data,
        std::vector<MorphingStrategy *> &strategies,
        ChemOperSelector *opers,
        RDKit::RWMol **newMols,
        std::string *smiles,
        std::string *formulas,
        double *weights,
        double *sascores,
        tbb::atomic<unsigned int> &kekulizeFailureCount,
        tbb::atomic<unsigned int> &sanitizeFailureCount,
        tbb::atomic<unsigned int> &morphingFailureCount
        );

    void operator()(const tbb::blocked_range<int> &r) const;

private:
    MorphingData &mData;
    std::vector<MorphingStrategy *> &mStrategies;

    ChemOperSelector *mOpers;
    RDKit::RWMol **mNewMols;
    std::string *mSmiles;
    std::string *mFormulas;
    double *mWeights;
    double *mSascore;
    tbb::atomic<unsigned int> &mKekulizeFailureCount;
    tbb::atomic<unsigned int> &mSanitizeFailureCount;
    tbb::atomic<unsigned int> &mMorphingFailureCount;
};
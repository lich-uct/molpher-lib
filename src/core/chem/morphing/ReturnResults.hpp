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

#include "chemoper_selectors.h"
#include "chem/fingerprintStrategy/FingerprintStrategy.h"
#include "chem/simCoefStrategy/SimCoefStrategy.h"
#include "MolpherMolecule.h"
#include "MorphingData.h"
#include "chem/morphingStrategy/MorphingStrategy.h"
#include "chem/SimCoefCalculator.hpp"

class ReturnResults
{
public:
    ReturnResults(
        RDKit::RWMol **newMols,
        std::string *smiles,
        std::string *formulas,
        std::string &parentSmile,
        ChemOperSelector *opers,
        double *weights,
        double *sascores,
        double *distToTarget,
        double *distToClosestDecoy,
        void *callerState,
        void (*deliver)(MolpherMolecule *, void *)
        );

    void operator()(const tbb::blocked_range<int> &r) const;

private:
    RDKit::RWMol **mNewMols;
    std::string *mSmiles;
    std::string *mFormulas;
    std::string &mParentSmile;
    ChemOperSelector *mOpers;
    double *mWeights;
    double *mSascore;
    double *mDistToTarget;
    double *mDistToClosestDecoy;

    void *mCallerState;
    void (*mDeliver)(MolpherMolecule *, void *);
    
    /**
     * Value of last next target decoy for new molecules. If the current molecule
     * distToClosestDecoy is equal to one then the value of lastDecoy for molecule
     * is increased by one up to mDecoySize.
     */
    int mNextDecoy;
    
    /**
     * The number of decoys used by the exploration process.
     */
    int mDecoySize;

};
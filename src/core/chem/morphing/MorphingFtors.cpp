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

#include <cfloat>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

//#include "main.hpp"
#include "core/misc/SynchRand.h"
#include "core/chem/ChemicalAuxiliary.h"
#include "core/chem/morphing/MorphingFtors.hpp"
#include "core/misc/SAScore.h"

/*
 Added code for working with SAScore
 */

CalculateMorphs::CalculateMorphs(
    MorphingData &data,
    std::vector<MorphingStrategy *> &strategies,
    ChemOperSelector *opers,
    RDKit::RWMol **newMols,
    std::string *smiles,
    std::string *formulas,
    double *weights,
    double *sascore, // added for SAScore
    tbb::atomic<unsigned int> &kekulizeFailureCount,
    tbb::atomic<unsigned int> &sanitizeFailureCount,
    tbb::atomic<unsigned int> &morphingFailureCount
    ) :
    mData(data),
    mStrategies(strategies),
    mOpers(opers),
    mNewMols(newMols),
    mSmiles(smiles),
    mFormulas(formulas),
    mWeights(weights),
    mSascore(sascore), // added for SAScore
    mKekulizeFailureCount(kekulizeFailureCount),
    mSanitizeFailureCount(sanitizeFailureCount),
    mMorphingFailureCount(morphingFailureCount)
{
    // no-op
}

void CalculateMorphs::operator()(const tbb::blocked_range<int> &r) const
{
//    DEBUG_REPORT("CalculateMorphs::operator")
    
    for (int i = r.begin(); i != r.end(); ++i) {
        int randPos = SynchRand::GetRandomNumber(mStrategies.size() - 1);
        MorphingStrategy *strategy = mStrategies[randPos];
        mOpers[i] = strategy->GetSelector();

        try {
            strategy->Morph(mData, &mNewMols[i]);
        } catch (const std::exception &exc) {
            ++mMorphingFailureCount; // atomic
            delete mNewMols[i];
            mNewMols[i] = NULL;
        }

        if (mNewMols[i]) {
            try {
                RDKit::RWMol *copy = new RDKit::RWMol();
                CopyMol(*(mNewMols[i]), *copy);
                delete mNewMols[i];
                mNewMols[i] = copy;

                mNewMols[i]->clearComputedProps();
                RDKit::MolOps::cleanUp(*(mNewMols[i]));
                mNewMols[i]->updatePropertyCache();
                RDKit::MolOps::Kekulize(*(mNewMols[i]));
                RDKit::MolOps::adjustHs(*(mNewMols[i]));

                //RDKit::MolOps::sanitizeMol(*(mNewMols[i]));

                mSmiles[i] = RDKit::MolToSmiles(*(mNewMols[i]));
                mFormulas[i] = RDKit::Descriptors::calcMolFormula(*(mNewMols[i]));
                mWeights[i] = RDKit::Descriptors::calcExactMW(*(mNewMols[i]));
                mSascore[i] = SAScore::getInstance()->getScore(*(mNewMols[i])); // added for SAScore
            } catch (const ValueErrorException &exc) {
                ++mKekulizeFailureCount; // atomic
                delete mNewMols[i];
                mNewMols[i] = NULL;
            } catch (const RDKit::MolSanitizeException &exc) {
                ++mSanitizeFailureCount; // atomic
                delete mNewMols[i];
                mNewMols[i] = NULL;
            } catch (const std::exception &exc) {
                ++mSanitizeFailureCount; // atomic
                delete mNewMols[i];
                mNewMols[i] = NULL;
            }
        }
    }
}

CalculateDistances::CalculateDistances(
    RDKit::RWMol **newMols,
    SimCoefCalculator &scCalc,
    Fingerprint *targetFp,
    std::vector<Fingerprint *> &decoysFp,
    double *distToTarget,
    double *distToClosestDecoy,
    int nextDecoy
    ) :
    mNewMols(newMols),
    mScCalc(scCalc),
    mTargetFp(targetFp),
    mDecoysFp(decoysFp),
    mDistToTarget(distToTarget),
    mDistToClosestDecoy(distToClosestDecoy),
    mNextDecoy(nextDecoy)
{
    // no-op
}

void CalculateDistances::operator()(const tbb::blocked_range<int> &r) const
{
    Fingerprint *fp;
    double dist, minDist;

    for (int i = r.begin(); i != r.end(); ++i) {
        if (mNewMols[i]) {
            fp = mScCalc.GetFingerprint(mNewMols[i]);
            mDistToTarget[i] = mScCalc.ConvertToDistance(
                mScCalc.GetSimCoef(mTargetFp, fp));

            dist = 0;
            // Calculate distance to the current decoy (mLastDecoy)
            if (mNextDecoy == -1 || mNextDecoy >= mDecoysFp.size() ) {
                // no calculation need, all the decoys are behind us 
                dist = 0;
            } else {
                // calculate distance to the next decoy for visit
                dist = mScCalc.ConvertToDistance(
                    mScCalc.GetSimCoef(mDecoysFp[mNextDecoy], fp));                
            }
            mDistToClosestDecoy[i] = dist;
            
            delete fp;
        }
    }
}

ReturnResults::ReturnResults(
    RDKit::RWMol **newMols,
    std::string *smiles,
    std::string *formulas,
    std::string &parentSmile,
    ChemOperSelector *opers,
    double *weights,
    double *sascore, // added for SAScore
    double *distToTarget,
    double *distToClosestDecoy,
    void *callerState,
    void (*deliver)(MolpherMolecule *, void *)
    ) :
    mNewMols(newMols),
    mSmiles(smiles),
    mFormulas(formulas),
    mParentSmile(parentSmile),
    mOpers(opers),
    mWeights(weights),
    mSascore(sascore), // added for SAScore
    mDistToTarget(distToTarget),
    mDistToClosestDecoy(distToClosestDecoy),
    mCallerState(callerState),
    mDeliver(deliver)
{
    // no-op
}

void ReturnResults::operator()(const tbb::blocked_range<int> &r) const
{
    // calculate value of the nextDecoy for molecules
    // that are close enough to the current decoy
    int nextDecoyForPassed = mNextDecoy;
    nextDecoyForPassed += mNextDecoy <= mDecoySize ? 1 : 0;
    
    for (int i = r.begin(); i != r.end(); ++i) {
        if (mNewMols[i]) {
            MolpherMolecule result(mSmiles[i], mFormulas[i], mParentSmile,
                mOpers[i], mDistToTarget[i], mDistToClosestDecoy[i],
                mWeights[i], mSascore[i]);
            
            /* Advance decoy functionality
            // are we close enough ?            
            if (result.distToClosestDecoy < mDecoyRange) {
                ++result.nextDecoy;
            }
            // update result.nextDecoy if need
            if (result.nextDecoy >= mDecoySize) {
                result.nextDecoy = mDecoySize;
            }*/

            mDeliver(&result, mCallerState);
        }
    }
}

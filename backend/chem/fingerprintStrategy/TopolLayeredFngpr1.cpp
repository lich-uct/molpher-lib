/* 
 * File:   TopolLayeredFngpr1.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "TopolLayeredFngpr1.hpp"

TopolLayeredFngpr1::TopolLayeredFngpr1(
    unsigned int layerFlags,
    unsigned int minPath,
    unsigned int maxPath,
    unsigned int fpSize,
    double tgtDensity,
    unsigned int minSize,
    std::vector<unsigned int> *atomCounts,
    ExplicitBitVect *setOnlyBits,
    bool branchedPaths
    ) :
    mLayerFlags(layerFlags),
    mMinPath(minPath),
    mMaxPath(maxPath),
    mFpSize(fpSize),
    mTgtDensity(tgtDensity),
    mMinSize(minSize),
    mAtomCounts(atomCounts),
    mSetOnlyBits(setOnlyBits),
    mBranchedPaths(branchedPaths)
{
     // no-op
}

TopolLayeredFngpr1::~TopolLayeredFngpr1()
{
    delete mAtomCounts;
    delete mSetOnlyBits;
}

Fingerprint *TopolLayeredFngpr1::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::LayeredFingerprintMol(*mol, mLayerFlags, mMinPath, mMaxPath,
        mFpSize, mTgtDensity, mMinSize, mAtomCounts, mSetOnlyBits,
        mBranchedPaths);
}

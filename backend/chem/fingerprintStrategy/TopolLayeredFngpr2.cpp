/* 
 * File:   TopolLayeredFngpr2.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "TopolLayeredFngpr2.hpp"

TopolLayeredFngpr2::TopolLayeredFngpr2(
    unsigned int layerFlags,
    unsigned int minPath,
    unsigned int maxPath,
    unsigned int fpSize,
    std::vector<unsigned int> *atomCounts,
    ExplicitBitVect *setOnlyBits,
    bool branchedPaths
    ) :
    mLayerFlags(layerFlags),
    mMinPath(minPath),
    mMaxPath(maxPath),
    mFpSize(fpSize),
    mAtomCounts(atomCounts),
    mSetOnlyBits(setOnlyBits),
    mBranchedPaths(branchedPaths)
{
     // no-op
}

TopolLayeredFngpr2::~TopolLayeredFngpr2()
{
    delete mAtomCounts;
    delete mSetOnlyBits;
}

Fingerprint *TopolLayeredFngpr2::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::LayeredFingerprintMol2(*mol, mLayerFlags, mMinPath, mMaxPath,
        mFpSize, mAtomCounts, mSetOnlyBits, mBranchedPaths);
}
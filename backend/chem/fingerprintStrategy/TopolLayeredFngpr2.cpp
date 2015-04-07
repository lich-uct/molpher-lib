/* 
 * File:   TopolLayeredFngpr2.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "TopolLayeredFngpr2.hpp"

TopolLayeredFngpr2::TopolLayeredFngpr2(
    unsigned int layerFlags,
    unsigned int fpSize,
    std::vector<unsigned int> *atomCounts,
    ExplicitBitVect *setOnlyBits
    ) :
    mLayerFlags(layerFlags),
    mFpSize(fpSize),
    mAtomCounts(atomCounts),
    mSetOnlyBits(setOnlyBits)
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
    return RDKit::PatternFingerprintMol(*mol, mFpSize, mAtomCounts, mSetOnlyBits);
}
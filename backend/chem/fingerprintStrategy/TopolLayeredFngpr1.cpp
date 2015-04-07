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
    unsigned int fpSize
    ) :
    mLayerFlags(layerFlags),
    mMinPath(minPath),
    mMaxPath(maxPath),
    mFpSize(fpSize)
{
     // no-op
}

TopolLayeredFngpr1::~TopolLayeredFngpr1()
{
}

Fingerprint *TopolLayeredFngpr1::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::LayeredFingerprintMol(*mol, mLayerFlags, mMinPath, mMaxPath, mFpSize);
}
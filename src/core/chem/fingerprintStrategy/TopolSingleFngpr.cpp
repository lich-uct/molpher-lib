/* 
 * File:   TopolSingleFngpr.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "TopolSingleFngpr.hpp"

TopolSingleFngpr::TopolSingleFngpr(
    unsigned int minPath,
    unsigned int maxPath,
    unsigned int fpSize,
    unsigned int nBitsPerHash,
    bool useHs,
    double tgtDensity,
    unsigned int minSize,
    bool branchedPaths
    ) :
    mMinPath(minPath),
    mMaxPath(maxPath),
    mFpSize(fpSize),
    mNBitsPerHash(nBitsPerHash),
    mUseHs(useHs),
    mTgtDensity(tgtDensity),
    mMinSize(minSize),
    mBranchedPaths(branchedPaths)
{
     // no-op
}

TopolSingleFngpr::~TopolSingleFngpr()
{
     // no-op
}


Fingerprint *TopolSingleFngpr::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::RDKFingerprintMol(*mol, mMinPath, mMaxPath, mFpSize,
        mNBitsPerHash, mUseHs, mTgtDensity, mMinSize, mBranchedPaths);
}
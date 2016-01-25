/* 
 * File:   TopolTorsFngpr.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "TopolTorsFngpr.hpp"

TopolTorsFngpr::TopolTorsFngpr(
    unsigned int nBits,
    unsigned int targetSize,
    const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    unsigned int nBitsPerEntry
    ) :
    mNBits(nBits),
    mTargetSize(targetSize),
    mFromAtoms(fromAtoms),
    mIgnoreAtoms(ignoreAtoms),
    mNBitsPerEntry(nBitsPerEntry)
{
     // no-op
}

TopolTorsFngpr::~TopolTorsFngpr()
{
    delete mFromAtoms;
    delete mIgnoreAtoms;
}

Fingerprint *TopolTorsFngpr::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(
        *mol, mNBits, mTargetSize, mFromAtoms, mIgnoreAtoms, 0, mNBitsPerEntry);
}
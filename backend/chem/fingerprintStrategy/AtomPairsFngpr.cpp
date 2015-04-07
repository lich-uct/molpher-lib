/* 
 * File:   AtomPairsFngpr.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "AtomPairsFngpr.hpp"

AtomPairsFngpr::AtomPairsFngpr(
    unsigned int nBits,
    unsigned int minLength,
    unsigned int maxLength,
    const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms
    ) :
    mNBits(nBits),
    mMinLength(minLength),
    mMaxLength(maxLength),
    mFromAtoms(fromAtoms),
    mIgnoreAtoms(ignoreAtoms)
{
     // no-op
}

AtomPairsFngpr::~AtomPairsFngpr()
{
    delete mFromAtoms;
    delete mIgnoreAtoms;
}

Fingerprint *AtomPairsFngpr::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(*mol,
        mNBits, mMinLength, mMaxLength, mFromAtoms, mIgnoreAtoms);
}

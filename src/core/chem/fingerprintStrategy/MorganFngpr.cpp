/* 
 * File:   MorganFngpr.cpp
 * Author: Petyr
 *
 * Created on 3. duben 2013
 */

#include "MorganFngpr.hpp"

MorganFngpr::MorganFngpr(
    unsigned int radius,
    unsigned int nBits,
    std::vector<boost::uint32_t> *invariants,
    const std::vector<boost::uint32_t> *fromAtoms,
    bool useChirality,
    bool useBondTypes,
    bool onlyNonzeroInvariants,
    RDKit::MorganFingerprints::BitInfoMap *atomsSettingBits
    ) :
    mRadius(radius),
    mNBits(nBits),
    mInvariants(invariants),
    mFromAtoms(fromAtoms),
    mUseChirality(useChirality),
    mUseBondTypes(useBondTypes),
    mOnlyNonzeroInvariants(onlyNonzeroInvariants),
    mAtomsSettingBits(atomsSettingBits)
{
     // no-op
}

MorganFngpr::~MorganFngpr()
{
    delete mInvariants;
    delete mFromAtoms;
    delete mAtomsSettingBits;
}

Fingerprint *MorganFngpr::GetFingerprint(RDKit::ROMol *mol)
{
    return RDKit::MorganFingerprints::getFingerprintAsBitVect(*mol, mRadius,
        mNBits, mInvariants, mFromAtoms, mUseChirality, mUseBondTypes, mOnlyNonzeroInvariants,
        mAtomsSettingBits);
}
/*
 Copyright (c) 2012 Vladimir Fiklik

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

#include "test_utils.h"

double GetTanimotoSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return TanimotoSimilarity(*fp1, *fp2);
}

double ConvertTanimotoToDistance(double coef)
{
    return 1 - coef;
}

double ComputeDistance(const MolpherMolecule &mol1, const MolpherMolecule &mol2)
{

    // TODO: Convert MolpherMolecule ro ROMol
    RDKit::ROMol *rMol1 = RDKit::SmilesToMol(mol1.smile);
    RDKit::ROMol *rMol2 = RDKit::SmilesToMol(mol2.smile);

    MorganFngpr mfpr;

    Fingerprint *fp1 = mfpr.GetFingerprint(rMol1);
    Fingerprint *fp2 = mfpr.GetFingerprint(rMol2);

    double simCoeff = GetTanimotoSimCoef(fp1, fp2);

    delete fp1;
    delete fp2;
    delete rMol1;
    delete rMol2;

    return ConvertTanimotoToDistance(simCoeff);

}

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
        mNBits, mInvariants, mFromAtoms, mUseBondTypes, mOnlyNonzeroInvariants,
        mAtomsSettingBits);
}
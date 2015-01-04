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

#pragma once

#include "global_types.h"
#include "MolpherMolecule.h"

#include <DataStructs/BitOps.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

double GetTanimotoSimCoef(Fingerprint *fp1, Fingerprint *fp2);
double ConvertTanimotoToDistance(double coef);
double ComputeDistance(const MolpherMolecule &mol1, const MolpherMolecule &mol2);

class MorganFngpr
{
public:
    /**
        Returns the Morgan fingerprint for a molecule.

        @param radius [in] the number of iterations to grow the fingerprint
        @param nBits [in] the number of bits in the final fingerprint
        @param invariants [in] optional pointer to a set of atom invariants to
            be used. By default ECFP-type invariants are used (calculated by
            getConnectivityInvariants())
        @param fromAtoms [in] if this is provided, only the atoms in the vector
            will be used as centers in the fingerprint
        @param useChirality [in] if set, additional information will be added to
            the fingerprint when chiral atoms are discovered. This will cause
            \verbatim C[C@H](F)Cl, C[C@@H](F)Cl, and CC(F)Cl \endverbatim to
            generate different fingerprints.
        @param useBondTypes [in] if set, bond types will be included as part of
            the hash for calculating bits
        @param onlyNonzeroInvariants [in] if set, bits will only be set from
            atoms that have a nonzero invariant.
        @param atomsSettingBits [in] if nonzero, this will be used to return
            information about the atoms that set each particular bit. The keys
            are the map are bit ids, the values are lists of (atomId, radius)
            pairs.
     */
    MorganFngpr(
        unsigned int radius = 1,
        unsigned int nBits = 2048,
        std::vector<boost::uint32_t> *invariants = 0,
        const std::vector<boost::uint32_t> *fromAtoms = 0,
        bool useChirality = false,
        bool useBondTypes = true,
        bool onlyNonzeroInvariants = false,
        RDKit::MorganFingerprints::BitInfoMap *atomsSettingBits = 0
        );

    ~MorganFngpr();

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

private:
    unsigned int mRadius;
    unsigned int mNBits;
    std::vector<boost::uint32_t> *mInvariants;
    const std::vector<boost::uint32_t> *mFromAtoms;
    bool mUseChirality;
    bool mUseBondTypes;
    bool mOnlyNonzeroInvariants;
    RDKit::MorganFingerprints::BitInfoMap *mAtomsSettingBits;
};
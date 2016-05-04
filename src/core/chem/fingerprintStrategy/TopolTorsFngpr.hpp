/*
 Copyright (c) 2012 Peter Szepe

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

#include "FingerprintStrategy.h"

class TopolTorsFngpr : public FingerprintStrategy
{
public:
    /**
        Returns a hashed topological-torsion fingerprint for a molecule.

        @param nBits [in] number of bits to include in the fingerprint
        @param targetSize [in] the number of atoms to include in the "torsions"
        @param fromAtoms [in] if provided, only torsions that start or end at
            the specified atoms will be included in the fingerprint
        @param ignoreAtoms [in] if provided, any torsions that include the
            specified atoms will not be included in the fingerprint
        @param nBitsPerEntry [in] number of bits to use in simulating counts
     */
    TopolTorsFngpr(
        unsigned int nBits = 2048,
        unsigned int targetSize = 4,
        const std::vector<boost::uint32_t> *fromAtoms = 0,
        const std::vector<boost::uint32_t> *ignoreAtoms = 0,
        unsigned int nBitsPerEntry = 4
        );

    ~TopolTorsFngpr();

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

private:
    unsigned int mNBits;
    unsigned int mTargetSize;
    const std::vector<boost::uint32_t> *mFromAtoms;
    const std::vector<boost::uint32_t> *mIgnoreAtoms;
    unsigned int mNBitsPerEntry;
};
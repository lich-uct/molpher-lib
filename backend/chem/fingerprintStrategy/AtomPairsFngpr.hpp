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

class AtomPairsFngpr : public FingerprintStrategy
{
public:
    /**
        Returns the hashed atom-pair fingerprint for a molecule.

        @param nBits [in] the length of the fingerprint to generate
        @param minLength [in] minimum distance between atoms to be  considered
            in a pair. Default is 1 bond.
        @param maxLength [in] maximum distance between atoms to be considered
            in a pair. Default is maxPathLen-1 bonds.
        @param fromAtoms [in] if provided, only atom pairs that involve the
            specified atoms will be included in the fingerprint
        @param ignoreAtoms [in] if provided, any atom pairs that include the
            specified atoms will not be included in the fingerprint
        @param nBitsPerEntry [in] number of bits to use in simulating counts
     */
    AtomPairsFngpr(
        unsigned int nBits = 2048,
        unsigned int minLength = 1,
        unsigned int maxLength = RDKit::AtomPairs::maxPathLen - 1,
        const std::vector<boost::uint32_t> *fromAtoms = 0,
        const std::vector<boost::uint32_t> *ignoreAtoms = 0,
        unsigned int nBitsPerEntry = 4
        );

    ~AtomPairsFngpr();

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

private:
    unsigned int mNBits;
    unsigned int mMinLength;
    unsigned int mMaxLength;
    const std::vector<boost::uint32_t> *mFromAtoms;
    const std::vector<boost::uint32_t> *mIgnoreAtoms;
    unsigned int mNBitsPerEntry;
};
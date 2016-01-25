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

class TopolSingleFngpr : public FingerprintStrategy
{
public:
    /**
        Generates a topological (Daylight like) fingerprint for a
        molecule using an alternate (faster) hashing algorithm

        @param minPath [in] the minimum path length (in bonds) to be included
        @param maxPath [in] the minimum path length (in bonds) to be included
        @param fpSize [in] the size of the fingerprint
        @param nBitsPerHash [in] the number of bits to be set by each path
        @param useHs [in] toggles inclusion of Hs in paths (if the molecule
            has explicit Hs)
        @param tgtDensity [in] if the generated fingerprint is below this
            density, it will be folded until the density is reached.
        @param minSize [in] the minimum size to which the fingerprint will be
            folded
        @param branchedPaths [in] toggles generation of branched subgraphs, not
            just linear paths
     */
    TopolSingleFngpr(
        unsigned int minPath = 1,
        unsigned int maxPath = 7,
        unsigned int fpSize = 2048,
        unsigned int nBitsPerHash = 2,
        bool useHs = true,
        double tgtDensity = 0.0,
        unsigned int minSize = 128,
        bool branchedPaths = true
        );

    ~TopolSingleFngpr();

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

private:
    unsigned int mMinPath;
    unsigned int mMaxPath;
    unsigned int mFpSize;
    unsigned int mNBitsPerHash;
    bool mUseHs;
    double mTgtDensity;
    unsigned int mMinSize;
    bool mBranchedPaths;
};
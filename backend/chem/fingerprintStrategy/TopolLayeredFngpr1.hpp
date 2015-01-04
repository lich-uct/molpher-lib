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

class TopolLayeredFngpr1 : public FingerprintStrategy
{
public:
    /**
        Generates a topological (Daylight like) fingerprint for a molecule using
        a layer-based hashing algorithm
        <b>Experimental:</b> This function is experimental. The API or results
        may change from release to release.

        @param layerFlags [in] the layers to be included (see below)
        @param minPath [in] the minimum path length (in bonds) to be included
        @param maxPath [in] the minimum path length (in bonds) to be included
        @param fpSize [in] the size of the fingerprint
        @param tgtDensity [in] if the generated fingerprint is below this
            density, it will be folded until the density is reached.
        @param minSize [in] the minimum size to which the fingerprint will be
            folded
        @param atomCounts [in] if provided, this will be used to provide the
            count of the number of paths that set bits each atom is involved in.
            The vector should have at least as many entries as the molecule has
            atoms and is not zeroed out here.
        @param setOnlyBits [in] if provided, only bits that are set in this bit
            vector will be set in the result. This is essentially the same as
            doing: (*res) &= (*setOnlyBits); but also has an impact on the
            atomCounts (if being used)
        @param branchedPaths [in] toggles generation of branched subgraphs, not
            just linear paths

        <b>Layer definitions:</b>
            - 0x01: pure topology
            - 0x02: bond order
            - 0x04: atom types
            - 0x08: presence of rings
            - 0x10: ring sizes
            - 0x20: aromaticity
     */
    TopolLayeredFngpr1(
        unsigned int layerFlags = 0xFFFFFFFF,
        unsigned int minPath = 1,
        unsigned int maxPath = 7,
        unsigned int fpSize = 2048,
        double tgtDensity = 0.0,
        unsigned int minSize = 128,
        std::vector<unsigned int> *atomCounts = 0,
        ExplicitBitVect *setOnlyBits = 0,
        bool branchedPaths = true
        );

    ~TopolLayeredFngpr1();

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

private:
    unsigned int mLayerFlags;
    unsigned int mMinPath;
    unsigned int mMaxPath;
    unsigned int mFpSize;
    double mTgtDensity;
    unsigned int mMinSize;
    std::vector<unsigned int> *mAtomCounts;
    ExplicitBitVect *mSetOnlyBits;
    bool mBranchedPaths;
};
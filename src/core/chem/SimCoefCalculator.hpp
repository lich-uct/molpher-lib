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

#include "core/chem/simCoefStrategy/SimCoefStrategy.h"

class SimCoefCalculator
{
    // the type of extended fingerprint building blocks
    typedef unsigned short extFpPart;

public:
    SimCoefCalculator(
        SimCoeffSelector sc = DEFAULT_SC,
        FingerprintSelector fp = DEFAULT_FP,
        RDKit::ROMol *source = NULL,
        RDKit::ROMol *target = NULL
        );
    ~SimCoefCalculator();

    double GetSimCoef(Fingerprint *fp1, Fingerprint *fp2);
    double GetSimCoef(RDKit::ROMol *mol1, RDKit::ROMol *mol2);
    double GetSimCoef(Fingerprint *fp1, RDKit::ROMol *mol2);

    double ConvertToDistance(double coef) const;

    Fingerprint *GetFingerprint(RDKit::ROMol *mol);

protected:
    Fingerprint *Extend(RDKit::ROMol *mol, Fingerprint *fp);

private:
    bool mExtended;
    std::map<AtomicNum, unsigned short> mAtomTypesToIdx;
    SimCoefStrategy *mScStrategy;
    FingerprintStrategy *mFpStrategy;
};
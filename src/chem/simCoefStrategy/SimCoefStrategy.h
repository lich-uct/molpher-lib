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

#include <map>

#include <DataStructs/BitOps.h>

#include "global_types.h"
#include "fingerprint_selectors.h"
#include "simcoeff_selectors.h"
#include "chem/fingerprintStrategy/FingerprintStrategy.h"

/*
    Legend:
    fp1_n: number of bits in vector 1
    fp1_o: number of on bits in vector 1
    (fp1&fp2)_o: number of on bits in the intersection of vectors 1 and 2
 */
class SimCoefStrategy
{
public:
    virtual double GetSimCoef(Fingerprint *fp1, Fingerprint *fp2) = 0;
    virtual double ConvertToDistance(double coef) = 0;
};

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

#include "SimCoefStrategy.h"

/// Tversky coefficients
/**
    (fp1&fp2)_o / [a*fp1_o + b*fp2_o + (1 - a - b)*(fp1&fp2)_o]  \n

    Notes: \n
    0 <= a, b <= 1 \n
    Tversky(a = 1, b = 1) = Tanimoto \n
    Tversky(a = 1 / 2, b = 1 / 2) = Dice \n
 */
class TverskySimCoef : public SimCoefStrategy
{
public:
    /// With no arguments, behaves same as Tanimoto.
    TverskySimCoef();

    /// Arguments are trimmed to belong into [0,1] interval.
    TverskySimCoef(double a, double b);

    double GetSimCoef(Fingerprint *fp1, Fingerprint *fp2);
    double ConvertToDistance(double coef);

private:
    double mA, mB;
};
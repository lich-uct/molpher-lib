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

/// Cosine coefficients
/**
    (fp1&fp2)_o / sqrt(fp1_o + fp2_o)
 */
class CosineSimCoef : public SimCoefStrategy
{
public:
    double GetSimCoef(Fingerprint *fp1, Fingerprint *fp2);
    double ConvertToDistance(double coef);
};
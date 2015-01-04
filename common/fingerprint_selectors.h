/*
 Copyright (c) 2012 Peter Szepe
 Copyright (c) 2012 Petr Koupy

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

#define DEFAULT_FP FP_MORGAN
#define MAX_STANDARD_FP FP_TOPOLOGICAL_TORSION

#include <string>

enum FingerprintSelector {
    FP_ATOM_PAIRS,
    FP_MORGAN,
    FP_TOPOLOGICAL,
    FP_TOPOLOGICAL_LAYERED_1,
    FP_TOPOLOGICAL_LAYERED_2,
    FP_TOPOLOGICAL_TORSION,
    FP_EXT_ATOM_PAIRS,
    FP_EXT_MORGAN,
    FP_EXT_TOPOLOGICAL,
    FP_EXT_TOPOLOGICAL_LAYERED_1,
    FP_EXT_TOPOLOGICAL_LAYERED_2,
    FP_EXT_TOPOLOGICAL_TORSION
};

const char *FingerprintShortDesc(const int selector);
const char *FingerprintLongDesc(const int selector);

FingerprintSelector FingerprintParse(const std::string& name);

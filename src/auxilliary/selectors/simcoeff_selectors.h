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

#define DEFAULT_SC SC_TANIMOTO

#include <string>

enum SimCoeffSelector {
    SC_ALL_BIT,
    SC_ASYMMETRIC,
    SC_BRAUN_BLANQUET,
    SC_COSINE,
    SC_DICE,
    SC_KULCZYNSKI,
    SC_MC_CONNAUGHEY,
    SC_ON_BIT,
    SC_RUSSEL,
    SC_SOKAL,
    SC_TANIMOTO,
    SC_TVERSKY_SUBSTRUCTURE,
    SC_TVERSKY_SUPERSTRUCTURE
};

const char *SimCoeffShortDesc(const int selector);
const char *SimCoeffLongDesc(const int selector);

SimCoeffSelector SimCoeffParse(const std::string& name);

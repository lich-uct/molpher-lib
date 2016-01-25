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

enum ChemOperSelector {
    OP_ADD_ATOM,
    OP_REMOVE_ATOM,
    OP_ADD_BOND,
    OP_REMOVE_BOND,
    OP_MUTATE_ATOM,
    OP_INTERLAY_ATOM,
    OP_BOND_REROUTE,
    OP_BOND_CONTRACTION
};

const char *ChemOperShortDesc(const int selector);
const char *ChemOperLongDesc(const int selector);

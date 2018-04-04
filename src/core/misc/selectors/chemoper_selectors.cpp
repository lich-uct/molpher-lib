/*
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

#include "selectors/chemoper_selectors.h"

static const char *shortDesc[] = {
    "OP_ADD_ATOM",
    "OP_REMOVE_ATOM",
    "OP_ADD_BOND",
    "OP_REMOVE_BOND",
    "OP_MUTATE_ATOM",
    "OP_INTERLAY_ATOM",
    "OP_BOND_REROUTE",
    "OP_BOND_CONTRACTION"
};

static const char *longDesc[] = {
    "Add Atom",
    "Remove Atom",
    "Add Bond",
    "Remove Bond",
    "Mutate Atom",
    "Interlay Atom",
    "Bond Reroute",
    "Bond Contraction"
};

std::string ChemOperShortDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(shortDesc) / sizeof(shortDesc[0])));
    if (validSelector) {
        return std::string(shortDesc[selector]);
    } else {
        return "";
    }
}

std::string ChemOperLongDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(longDesc) / sizeof(longDesc[0])));
    if (validSelector) {
        return std::string(longDesc[selector]);
    } else {
        return "";
    }
}

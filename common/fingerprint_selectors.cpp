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

#include "fingerprint_selectors.h"

#include <stdexcept>
#include "boost/algorithm/string/predicate.hpp"

static const char *shortDesc[] = {
    "ATP",
    "MRG",
    "TOP",
    "TL1",
    "TL2",
    "TPT",
    "EATP",
    "EMRG",
    "ETOP",
    "ETL1",
    "ETL2",
    "ETPT"
};

static const char *longDesc[] = {
    "Atom Pairs",
    "Morgan",
    "Topological",
    "Topological Layered 1",
    "Topological Layered 2",
    "Topological Torsion",
    "Atom Pairs (ext)",
    "Morgan (ext)",
    "Topological (ext)",
    "Topological Layered 1 (ext)",
    "Topological Layered 2 (ext)",
    "Topological Torsion (ext)"
};

const char *FingerprintShortDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(shortDesc) / sizeof(shortDesc[0])));
    if (validSelector) {
        return shortDesc[selector];
    } else {
        return "";
    }
}

const char *FingerprintLongDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(longDesc) / sizeof(longDesc[0])));
    if (validSelector) {
        return longDesc[selector];
    } else {
        return "";
    }
}

FingerprintSelector FingerprintParse(const std::string& name) {
    if (boost::iequals(name, "FP_ATOM_PAIRS")) {
        return FP_ATOM_PAIRS;
    } if (boost::iequals(name, "FP_MORGAN")) {
        return FP_MORGAN;
    } if (boost::iequals(name, "FP_TOPOLOGICAL")) {
        return FP_TOPOLOGICAL;
    } if (boost::iequals(name, "FP_TOPOLOGICAL_LAYERED_1")) {
        return FP_TOPOLOGICAL_LAYERED_1;
    } if (boost::iequals(name, "FP_TOPOLOGICAL_LAYERED_2")) {
        return FP_TOPOLOGICAL_LAYERED_2;
    } if (boost::iequals(name, "FP_TOPOLOGICAL_TORSION")) {
        return FP_TOPOLOGICAL_TORSION;
    } if (boost::iequals(name, "FP_EXT_ATOM_PAIRS")) {
        return FP_EXT_ATOM_PAIRS;
    } if (boost::iequals(name, "FP_EXT_MORGAN")) {
        return FP_EXT_MORGAN;
    } if (boost::iequals(name, "FP_EXT_TOPOLOGICAL")) {
        return FP_EXT_TOPOLOGICAL;
    } if (boost::iequals(name, "FP_EXT_TOPOLOGICAL_LAYERED_1")) {
        return FP_EXT_TOPOLOGICAL_LAYERED_1;
    } if (boost::iequals(name, "FP_EXT_TOPOLOGICAL_LAYERED_2")) {
        return FP_EXT_TOPOLOGICAL_LAYERED_2;
    } if (boost::iequals(name, "FP_EXT_TOPOLOGICAL_TORSION")) {
        return FP_EXT_TOPOLOGICAL_TORSION;
    } else {
        throw std::runtime_error("Unknown fingerprint name.");
    }    
}
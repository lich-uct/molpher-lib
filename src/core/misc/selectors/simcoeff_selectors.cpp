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

#include "simcoeff_selectors.h"

#include <stdexcept>
#include "boost/algorithm/string/predicate.hpp"

static const char *shortDesc[] = {
    "ALB",
    "ASM",
    "BBL",
    "COS",
    "DIC",
    "KLC",
    "MCC",
    "ONB",
    "RSL",
    "SKL",
    "TNM",
    "TSB",
    "TSP"
};

static const char *longDesc[] = {
    "All Bit",
    "Asymmetric",
    "Braun-Blanquet",
    "Cosine",
    "Dice",
    "Kulczynski",
    "McConnaughey",
    "On Bit",
    "Russel",
    "Sokal",
    "Tanimoto",
    "Tversky (substructure)",
    "Tversky (superstructure)"
};

const char *SimCoeffShortDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(shortDesc) / sizeof(shortDesc[0])));
    if (validSelector) {
        return shortDesc[selector];
    } else {
        return "";
    }
}

const char *SimCoeffLongDesc(const int selector)
{
    bool validSelector = (selector >= 0) &&
        (selector < (int)(sizeof(longDesc) / sizeof(longDesc[0])));
    if (validSelector) {
        return longDesc[selector];
    } else {
        return "";
    }
}

SimCoeffSelector SimCoeffParse(const std::string& name) {
    if (boost::iequals(name, "SC_ALL_BIT")) {
        return SC_ALL_BIT;
    } if (boost::iequals(name, "SC_ASYMMETRIC")) {
        return SC_ASYMMETRIC;
    } if (boost::iequals(name, "SC_BRAUN_BLANQUET")) {
        return SC_BRAUN_BLANQUET;
    } if (boost::iequals(name, "SC_COSINE")) {
        return SC_COSINE;
    } if (boost::iequals(name, "SC_DICE")) {
        return SC_DICE;
    } if (boost::iequals(name, "SC_KULCZYNSKI")) {
        return SC_KULCZYNSKI;
    } if (boost::iequals(name, "SC_MC_CONNAUGHEY")) {
        return SC_MC_CONNAUGHEY;
    } if (boost::iequals(name, "SC_ON_BIT")) {
        return SC_ON_BIT;
    } if (boost::iequals(name, "SC_RUSSEL")) {
        return SC_RUSSEL;
    } if (boost::iequals(name, "SC_SOKAL")) {
        return SC_SOKAL;
    } if (boost::iequals(name, "SC_TANIMOTO")) {
        return SC_TANIMOTO;
    } if (boost::iequals(name, "SC_TVERSKY_SUBSTRUCTURE")) {
        return SC_TVERSKY_SUBSTRUCTURE;
        } if (boost::iequals(name, "SC_TVERSKY_SUPERSTRUCTURE")) {
        return SC_TVERSKY_SUPERSTRUCTURE;
    } else {
        throw std::runtime_error("Unknown fingerprint name.");
    }    
}

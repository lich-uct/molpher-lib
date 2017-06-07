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

#include <GraphMol/GraphMol.h>

#include "core/misc/global_types.h"
#include "selectors/chemoper_selectors.h"
#include "core/chem/morphing/MorphingData.h"

class MorphingStrategy
{
public:
    virtual void Morph(MorphingData &data, RDKit::RWMol **nMol) = 0;
    virtual ChemOperSelector GetSelector() = 0;
};

void InitStrategies(
		const std::vector<ChemOperSelector> &chemOperSelectors,
		std::vector<MorphingStrategy *> &strategies);

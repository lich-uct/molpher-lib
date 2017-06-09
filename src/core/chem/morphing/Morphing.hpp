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

#include <vector>

#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"
#include "selectors/chemoper_selectors.h"
#include "data_structs/MolpherMol.hpp"

#ifndef MORPHING_REPORTING
#define MORPHING_REPORTING 1
#endif

// TODO: get rid of this function and handle all morphing through MolpherMol instances
void GenerateMorphs(
    MolpherMol &candidate,
    unsigned int morphAttempts,
    FingerprintSelector fingerprintSelector,
    SimCoeffSelector simCoeffSelector,
    std::vector<ChemOperSelector> &chemOperSelectors,
    MolpherMol &target,
    std::vector<MolpherMol> &decoys,
    tbb::task_group_context &tbbCtx,
    void *callerState,
    void (*deliver)(std::shared_ptr<MolpherMol>, void *)
    );

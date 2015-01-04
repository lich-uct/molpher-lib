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

#include "PathFinderContext.h"

void PathFinderContext::ContextToSnapshot(
    const PathFinderContext &ctx, IterationSnapshot &snp)
{
    snp.jobId = ctx.jobId;
    snp.iterIdx = ctx.iterIdx;
    snp.elapsedSeconds = ctx.elapsedSeconds;

    snp.fingerprintSelector = ctx.fingerprintSelector;
    snp.simCoeffSelector = ctx.simCoeffSelector;
    snp.dimRedSelector = ctx.dimRedSelector;

    snp.chemOperSelectors.clear();
    snp.chemOperSelectors.resize(ctx.chemOperSelectors.size(), 0);
    for (size_t i = 0; i < ctx.chemOperSelectors.size(); ++i) {
        snp.chemOperSelectors[i] = ctx.chemOperSelectors[i];
    }

    snp.params = ctx.params;

    snp.source = ctx.source;
    snp.target = ctx.target;
    snp.decoys = ctx.decoys;

    snp.candidates.clear();
    for (CandidateMap::const_iterator it = ctx.candidates.begin();
            it != ctx.candidates.end(); it++) {
        snp.candidates.insert(*it);
    }

    snp.morphDerivations.clear();
    for (MorphDerivationMap::const_iterator it = ctx.morphDerivations.begin();
            it != ctx.morphDerivations.end(); it++) {
        snp.morphDerivations.insert(*it);
    }

    snp.prunedDuringThisIter.clear();
    for (PrunedMoleculeVector::const_iterator it = ctx.prunedDuringThisIter.begin();
            it != ctx.prunedDuringThisIter.end(); it++) {
        snp.prunedDuringThisIter.push_back(*it);
    }
}

void PathFinderContext::SnapshotToContext(
    const IterationSnapshot &snp, PathFinderContext &ctx)
{
    ctx.jobId = snp.jobId;
    ctx.iterIdx = snp.iterIdx;
    ctx.elapsedSeconds = snp.elapsedSeconds;

    ctx.fingerprintSelector = (FingerprintSelector) snp.fingerprintSelector;
    ctx.simCoeffSelector = (SimCoeffSelector) snp.simCoeffSelector;
    ctx.dimRedSelector = (DimRedSelector) snp.dimRedSelector;

    ctx.chemOperSelectors.clear();
    ctx.chemOperSelectors.resize(snp.chemOperSelectors.size(), (ChemOperSelector) 0);
    for (size_t i = 0; i < snp.chemOperSelectors.size(); ++i) {
        ctx.chemOperSelectors[i] = (ChemOperSelector) snp.chemOperSelectors[i];
    }

    ctx.params = snp.params;

    ctx.source = snp.source;
    ctx.target = snp.target;
    ctx.decoys = snp.decoys;

    ctx.candidates.clear();
    for (IterationSnapshot::CandidateMap::const_iterator it = snp.candidates.begin();
            it != snp.candidates.end(); it++) {
        ctx.candidates.insert(*it);
    }

    ctx.morphDerivations.clear();
    for (IterationSnapshot::MorphDerivationMap::const_iterator it = snp.morphDerivations.begin();
            it != snp.morphDerivations.end(); it++) {
        ctx.morphDerivations.insert(*it);
    }

    ctx.prunedDuringThisIter.clear();
    ctx.prunedDuringThisIter.reserve(snp.prunedDuringThisIter.size());
    for (IterationSnapshot::PrunedMoleculeVector::const_iterator it = snp.prunedDuringThisIter.begin();
            it != snp.prunedDuringThisIter.end(); it++) {
        ctx.prunedDuringThisIter.push_back(*it);
    }
}

void PathFinderContext::ContextToLightSnapshot(
    const PathFinderContext &ctx, IterationSnapshot &snp)
{
    snp.jobId = ctx.jobId;
    snp.iterIdx = ctx.iterIdx;
    snp.elapsedSeconds = ctx.elapsedSeconds;

    snp.fingerprintSelector = ctx.fingerprintSelector;
    snp.simCoeffSelector = ctx.simCoeffSelector;
    snp.dimRedSelector = ctx.dimRedSelector;

    snp.chemOperSelectors.clear();
    snp.chemOperSelectors.resize(ctx.chemOperSelectors.size(), 0);
    for (size_t i = 0; i < ctx.chemOperSelectors.size(); ++i) {
        snp.chemOperSelectors[i] = ctx.chemOperSelectors[i];
    }

    snp.params = ctx.params;

    snp.source = ctx.source;
    snp.target = ctx.target;
    snp.decoys = ctx.decoys;
}

void PathFinderContext::clear()
{
    chemOperSelectors.clear();
    decoys.clear();
    candidates.clear();
    morphDerivations.clear();
    prunedDuringThisIter.clear();
    prunedDuringThisIter.shrink_to_fit();
}

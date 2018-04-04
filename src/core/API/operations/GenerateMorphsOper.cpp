/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2016 Martin Šícho

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

#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#include <morphing/operators/MorphingOperator.hpp>
#include <morphing/Molpher.hpp>
#include <core/misc/SAScore.h>
#include <core/misc/inout.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "operations/GenerateMorphsOper.hpp"
#include "GenerateMorphsOperImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"

GenerateMorphsOper::GenerateMorphsOper(std::shared_ptr<ExplorationTree> expTree, bool set_tree_ownership) : 
pimpl(new GenerateMorphsOper::GenerateMorphsOperImpl(expTree, set_tree_ownership)) 
{
    setTreeOperPimpl(pimpl);
}

GenerateMorphsOper::GenerateMorphsOper(bool set_tree_ownership) : 
pimpl(new GenerateMorphsOper::GenerateMorphsOperImpl(set_tree_ownership))
{
    setTreeOperPimpl(pimpl);
}

void GenerateMorphsOper::operator()() {
    (*pimpl)();
}

GenerateMorphsOper::GenerateMorphsOperImpl::GenerateMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree, bool set_tree_ownership) :
TreeOperation::TreeOperationImpl::TreeOperationImpl(expTree)
, mSetTreeOwnershipForMorphs(set_tree_ownership)
{
    // no action
}

GenerateMorphsOper::GenerateMorphsOperImpl::GenerateMorphsOperImpl(bool set_tree_ownership) : 
TreeOperation::TreeOperationImpl::TreeOperationImpl()
, mSetTreeOwnershipForMorphs(set_tree_ownership)
{
    // no action
}

//void CollectMorphs::MorphCollector(std::shared_ptr<MolpherMol> morph, void *functor) {
//    CollectMorphs *collect =
//            (CollectMorphs *) functor;
//    (*collect)(morph);
//}

CollectMorphs::CollectMorphs(
	ConcurrentMolVector &morphs
	, std::shared_ptr<ExplorationTree> tree
	, bool set_ownership
	, std::shared_ptr<MolpherMol> target
	, FingerprintSelector fing_sel
	, SimCoeffSelector sim_sel)
:
MorphCollector()
	, mMorphs(morphs)
	, mTree(tree)
	, mSetTreeOwnership(set_ownership)
	, target(target)
	, sim_select(sim_sel)
	, fing_select(fing_sel)
	, mScCalc(sim_sel, fing_sel)
	, mTargetFp(nullptr)
{
	mCollectAttemptCount = 0;
	if (target) {
		auto target_rd = target->asRDMol();
		mTargetFp.reset(mScCalc.GetFingerprint(target_rd));
		delete target_rd;
	}
}

CollectMorphs::CollectMorphs(
    ConcurrentMolVector &morphs
    , std::shared_ptr<ExplorationTree> tree
    , bool set_ownership) 
:
CollectMorphs(
	morphs
	, tree
	, set_ownership
	, nullptr
	, DEFAULT_FP
	, DEFAULT_SC
)
{
    // no action
}

CollectMorphs::CollectMorphs(ConcurrentMolVector &morphs)
:
CollectMorphs(
	morphs
	, nullptr
	, false
	, nullptr
	, DEFAULT_FP
	, DEFAULT_SC
)
{
	// no action
}

CollectMorphs::CollectMorphs(
		ConcurrentMolVector &morphs
		, std::shared_ptr<MolpherMol> target
		, FingerprintSelector fing_sel
		, SimCoeffSelector sim_sel)
:
CollectMorphs(
	morphs
	, nullptr
	, false
	, target
	, fing_sel
	, sim_sel
)
{
	// no action
}

void CollectMorphs::operator()(
		std::shared_ptr<MolpherMol> morph
		, std::shared_ptr<MorphingOperator> operator_
) {
	if (!morph) return; // ignore empty morph
    ++mCollectAttemptCount; // atomic
    ConcurrentSmileSet::const_accessor dummy;
    if (mDuplicateChecker.insert(dummy, morph->getSMILES())) {
        morph->setParentSMILES(operator_->getOriginal()->getSMILES());
		morph->setParentOper(operator_->getName());

		Fingerprint* fp(nullptr);
		RDKit::RWMol* morph_rd(nullptr);
		if (target) {
			try {
				morph_rd = morph->asRDMol();
				fp = mScCalc.GetFingerprint(morph_rd);
				morph->setDistToTarget(mScCalc.ConvertToDistance(
						mScCalc.GetSimCoef(mTargetFp.get(), fp)));
				delete morph_rd;
				delete fp;
			} catch (const std::exception &exc) {
				SynchCerr(
						"Failed to calculate fingerprint due to: " + std::string(exc.what())
						+ "\n\tParent molecule: " + operator_->getOriginal()->getSMILES()
						+ "\n\tParent operator: " + operator_->getName()
						+ "\n\tGenerated morph: " + morph->getSMILES()
				);
				if (fp) delete fp;
				if (morph_rd) delete morph_rd;
				return;
			}
		}
		if (mTree && mSetTreeOwnership) {
            morph->setOwner(mTree);
        }
        mMorphs.push_back(morph);
    } else {
        // ignore duplicate
    }
}

unsigned int CollectMorphs::WithdrawCollectAttemptCount() {
    unsigned int ret = mCollectAttemptCount;
    mCollectAttemptCount = 0;
    return ret;
}

void GenerateMorphsOper::GenerateMorphsOperImpl::operator()() {
    auto tree = getTree();
    if (tree) {
        auto tree_pimpl = tree->pimpl;

        ConcurrentMolVector leaves;
        tree_pimpl->fetchLeaves(tree, true, leaves);
        
        tree_pimpl->candidates.clear();
        std::shared_ptr<CollectMorphs> collector(new CollectMorphs(
				tree_pimpl->candidates
				, tree
				, mSetTreeOwnershipForMorphs
				, tree_pimpl->target
				, (FingerprintSelector) tree_pimpl->fingerprint
				, (SimCoeffSelector) tree_pimpl->simCoeff
		));
        for (auto& leaf : leaves) {
            unsigned int morphAttempts = tree_pimpl->params.cntMorphs;
            if (leaf->getDistToTarget() < tree_pimpl->params.distToTargetDepthSwitch) {
                morphAttempts = tree_pimpl->params.cntMorphsInDepth;
            }

			// generate morphs for the current leaf
			try {
				Molpher molpher(
						leaf
						, tree_pimpl->chemOpers
						, tree_pimpl->threadCnt
						, morphAttempts
						, std::static_pointer_cast<MorphCollector>(collector)
				);
				tree_pimpl->candidates.reserve(tree_pimpl->candidates.size() + morphAttempts);
				molpher();
			} catch (std::exception exp) {
				Cerr("Morphing failed for leaf: " + leaf->getSMILES());
				tree->deleteSubtree(leaf->getSMILES(), false);
				continue;
			}

            MorphDerivationMap::accessor ac;
            if (tree_pimpl->morphDerivations.find(ac, leaf->getSMILES())) {
                ac->second += collector->WithdrawCollectAttemptCount();
            } else {
                tree_pimpl->morphDerivations.insert(ac, leaf->getSMILES());
                ac->second = collector->WithdrawCollectAttemptCount();
            }
        }
        tree_pimpl->candidates.shrink_to_fit();
        tree_pimpl->candidatesMask.clear();
        tree_pimpl->candidatesMask.resize(tree_pimpl->candidates.size(), true);
//        SynchCout("Generated " + parseNumber(tree_pimpl->candidates.size()) + " morphs.");
    } else {
        throw std::runtime_error("Cannot generate morphs. No tree associated with this instance.");
    }
}

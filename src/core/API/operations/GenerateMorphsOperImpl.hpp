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

#ifndef GENERATEMORPHSOPERIMPL_HPP
#define	GENERATEMORPHSOPERIMPL_HPP

#include "core/misc/global_types.h"
#include "TreeOperationImpl.hpp"
#include "operations/GenerateMorphsOper.hpp"

class CollectMorphs
{
public:
    CollectMorphs(ConcurrentMolVector &morphs, std::shared_ptr<ExplorationTree> tree, bool set_ownership);
    CollectMorphs(ConcurrentMolVector &morphs, bool set_ownership); // TODO: remove the set_ownership parameter
    void operator()(std::shared_ptr<MolpherMol> morph);
    unsigned int WithdrawCollectAttemptCount();
    static void MorphCollector(std::shared_ptr<MolpherMol> morph, void *functor);

private:
    ConcurrentSmileSet mDuplicateChecker;
    ConcurrentMolVector &mMorphs;
    std::shared_ptr<ExplorationTree> mTree;
    bool mSetTreeOwnership;
    tbb::atomic<unsigned int> mCollectAttemptCount;
};

class GenerateMorphsOper::GenerateMorphsOperImpl : public TreeOperation::TreeOperationImpl {

private:
    bool mSetTreeOwnershipForMorphs;
    
public:
    GenerateMorphsOperImpl(std::shared_ptr<ExplorationTree> expTree, bool set_ownership = false);
    GenerateMorphsOperImpl(bool set_ownership = false);
    virtual void operator()();

};

#endif	/* GENERATEMORPHSOPERIMPL_HPP */


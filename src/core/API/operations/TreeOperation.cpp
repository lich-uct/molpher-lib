/*
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

#include "operations/TreeOperation.hpp"
#include "TreeOperationImpl.hpp"
#include "core/API/ExplorationTreeImpl.h"


TreeOperation::TreeOperation(std::shared_ptr<ExplorationTree> expTree) : 
pimpl(new TreeOperation::TreeOperationImpl(expTree))
{
    // no action
}

TreeOperation::TreeOperation() :
pimpl(new TreeOperation::TreeOperationImpl())
{
    // no action
}

TreeOperation::~TreeOperation() {
    // no action
}

std::shared_ptr<ExplorationTree> TreeOperation::getTree() {
    return pimpl->getTree();
}

void TreeOperation::setTree(std::shared_ptr<ExplorationTree> tree) {
    pimpl->setTree(tree);
}

TreeOperation::TreeOperationImpl::TreeOperationImpl(std::shared_ptr<ExplorationTree> expTree) :
tree(expTree)
{
    // no action
}

TreeOperation::TreeOperationImpl::TreeOperationImpl() {
    // no action
}

std::shared_ptr<ExplorationTree> TreeOperation::TreeOperationImpl::getTree() {
    return tree;
}

void TreeOperation::TreeOperationImpl::setTree(std::shared_ptr<ExplorationTree> tree) {
    this->tree = tree;
}

void TreeOperation::setTreeOperPimpl(std::shared_ptr<TreeOperationImpl> pimpl) {
    this->pimpl = pimpl;
}


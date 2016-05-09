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

%{    
// data structs
#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationData.hpp"
#include "data_structs/ExplorationTree.hpp"

// operations
#include "operations/TreeOperation.hpp"
#include "operations/FindLeavesOper.hpp"
#include "operations/GenerateMorphsOper.hpp"
#include "operations/SortMorphsOper.hpp"
#include "operations/FilterMorphsOper.hpp"
#include "operations/ExtendTreeOper.hpp"
#include "operations/PruneTreeOper.hpp"
#include "operations/TraverseOper.hpp"

// callbacks
#include "operations/callbacks/TraverseCallback.hpp"

// selectors
#include "selectors/chemoper_selectors.h"
#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"
%}

%template(StringSet) std::set<std::string>;
%template(IntSet) std::set<int>;
%template(StringVector) std::vector<std::string>;
%template(BoolVector) std::vector<bool>;
%template(MolpherMolVector) std::vector<std::shared_ptr<MolpherMol> >;
%template(MolpherMolMap) std::map<std::string, std::shared_ptr<MolpherMol> >;

%shared_ptr(ExplorationTree);
%shared_ptr(MolpherMol);
%shared_ptr(ExplorationData);
%shared_ptr(std::shared_ptr<std::map<std::string, std::shared_ptr<MolpherMol> > >);
%shared_ptr(std::shared_ptr<std::vector<std::shared_ptr<MolpherMol> > >);

// data structs
%include "data_structs/data_structs.i"

// tree operations
%include "operations/operations.i"

// callbacks
%include "operations/callbacks/callbacks.i"
        
// selectors
%include "selectors/selectors.i"
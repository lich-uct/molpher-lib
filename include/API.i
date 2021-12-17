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
#include "data_structs/MolpherAtom.hpp"
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
#include "operations/CleanMorphsOper.hpp"

// callbacks
#include "operations/callbacks/TraverseCallback.hpp"
#include "operations/callbacks/SortMorphsCallback.hpp"

// selectors
#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"

// morphing facilities
#include "morphing/operators/MorphingOperator.hpp"
#include "morphing/operators/AddAtom.hpp"
#include "morphing/operators/RemoveAtom.hpp"
#include "morphing/operators/MutateAtom.hpp"
#include "morphing/operators/AddBond.hpp"
#include "morphing/operators/RemoveBond.hpp"
#include "morphing/operators/ContractBond.hpp"
#include "morphing/operators/InterlayAtom.hpp"
#include "morphing/operators/RerouteBond.hpp"
#include "morphing/operators/ReactionOperator.hpp"
#include "morphing/AtomLibrary.hpp"
#include "morphing/Molpher.hpp"
#include "morphing/MorphCollector.hpp"
%}

%template(StringSet) std::set<std::string>;
%template(IntSet) std::set<int>;
%template(StringVector) std::vector<std::string>;
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(UIntVector) std::vector<unsigned int>;
%template(UIntVectorPairVector) std::vector<std::pair<unsigned int, unsigned int> >;
%template(BoolVector) std::vector<bool>;
%template(MolpherMolVector) std::vector<std::shared_ptr<MolpherMol> >;
%template(MolpherMolMap) std::map<std::string, std::shared_ptr<MolpherMol> >;
%template(MolpherAtomVector) std::vector<std::shared_ptr<MolpherAtom> >;
%template(MorphingOperatorVector) std::vector<std::shared_ptr<MorphingOperator> >;
%template(MorphCollectorVector) std::vector<std::shared_ptr<MorphCollector> >;

%shared_ptr(ExplorationTree);
%shared_ptr(MolpherMol);
%shared_ptr(MolpherAtom);
%shared_ptr(ExplorationData);
%shared_ptr(MorphingOperator);
%shared_ptr(AddAtom);
%shared_ptr(RemoveAtom);
%shared_ptr(MutateAtom);
%shared_ptr(AddBond);
%shared_ptr(RemoveBond);
%shared_ptr(ContractBond);
%shared_ptr(InterlayAtom);
%shared_ptr(RerouteBond);
%shared_ptr(ReactionOperator);
%shared_ptr(MorphCollector);
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

// morphing
%include "morphing/morphing.i"

// morphing operators
%include "morphing/operators/operators.i"
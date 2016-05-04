
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
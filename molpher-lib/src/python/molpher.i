%module molpher
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
#include "../../include/morphing_functions.hpp"
#include "../../include/molpher_API/MolpherMol.hpp"
#include "../../include/molpher_API/ExplorationParameters.hpp"
#include "../../include/molpher_API/ExplorationTreeSnapshot.hpp"
#include "../../include/molpher_API/ExplorationTree.hpp"
#include "extensions/SAScore.h"
#include "../../include/molpher_API/operations/TreeOperation.hpp"
#include "../../include/molpher_API/operations/FindLeavesOper.hpp"
#include "../../include/molpher_API/operations/GenerateMorphsOper.hpp"
#include "../../include/molpher_API/operations/SortMorphsOper.hpp"
#include "../../include/molpher_API/operations/FilterMorphsOper.hpp"
#include "../../include/molpher_API/operations/ExtendTreeOper.hpp"
%}

%init %{
    SAScore::loadData();
%}

// complete morphing function
%include "../../include/morphing_functions.hpp"

// tree operations
%include "../../include/molpher_API/operations/TreeOperation.hpp"
%include "../../include/molpher_API/operations/FindLeavesOper.hpp"
%include "../../include/molpher_API/operations/GenerateMorphsOper.hpp"
%include "../../include/molpher_API/operations/SortMorphsOper.hpp"
%include "../../include/molpher_API/operations/FilterMorphsOper.hpp"
%include "../../include/molpher_API/operations/ExtendTreeOper.hpp"

// MoplherMolecule wrapper
%ignore MolpherMol::MolpherMol(const MolpherMolecule& mol);
%ignore MolpherMol::getMol();
%ignore MolpherMol::MolpherMol();
%include "../../include/molpher_API/MolpherMol.hpp"

// ExplorationParameters wrapper
%include "../../include/molpher_API/ExplorationParameters.hpp"
        
// ExplorationTreeSnapshot wrapper
%include "../../include/molpher_API/ExplorationTreeSnapshot.hpp"
        
// ExplorationTree wrapper
%template(MolpherMolVector) std::vector<MolpherMol>;
%template(BoolVector) std::vector<bool>;
%catches(std::runtime_error) ExplorationTree::setCandidateMorphsMask(std::vector<bool>);
%ignore ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret);
%include "../../include/molpher_API/ExplorationTree.hpp";
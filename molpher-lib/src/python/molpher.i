%module molpher
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
#include "../../include/morphing_functions.hpp"
#include "../../include/molpher_API/MolpherMol.hpp"
#include "../../include/molpher_API/ExplorationParameters.hpp"
#include "../../include/molpher_API/ExplorationTreeSnapshot.hpp"
    #include "../../include/molpher_API/ExplorationTree.hpp"
%}



// complete morphing function
%include "../../include/morphing_functions.hpp"

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
%ignore ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret);
%include "../../include/molpher_API/ExplorationTree.hpp"
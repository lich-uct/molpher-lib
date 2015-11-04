
%{
#include "molpher_API/MolpherMol.hpp"
#include "molpher_API/ExplorationParameters.hpp"
#include "molpher_API/ExplorationTreeSnapshot.hpp"
#include "molpher_API/ExplorationTree.hpp"
%}

// tree operations
%include "operations/operations.i"

// callbacks
%include "callbacks/callbacks.i"

// the MoplherMolecule wrapper class
%ignore MolpherMol::MolpherMol(MolpherMolecule& mol);
%ignore MolpherMol::getMol();
%ignore MolpherMol::MolpherMol();
%ignore MolpherMol::operator=(const MolpherMol&);
%include "MolpherMol.hpp"

// ExplorationParameters wrapper
%include "ExplorationParameters.hpp"
        
// ExplorationTreeSnapshot wrapper
%newobject ExplorationTreeSnapshot::load(const std::string& filename);
%catches(std::runtime_error) ExplorationTreeSnapshot::load(const std::string& filename);
%include "ExplorationTreeSnapshot.hpp"
        
// ExplorationTree wrapper
%template(MolpherMolVector) std::vector<MolpherMol>;
%template(BoolVector) std::vector<bool>;
%newobject ExplorationTree::createSnapshot();
%newobject ExplorationTree::createFromSnapshot(ExplorationTreeSnapshot& snapshot);
//%newobject ExplorationTree::fetchLeaves();
%newobject ExplorationTree::fetchMol(const std::string& canonSMILES);
%newobject ExplorationTree::getCandidateMorphs();
%newobject ExplorationTree::getCandidateMorphsMask();
%catches(std::runtime_error) ExplorationTree::setCandidateMorphsMask(const std::vector<bool>&);
%catches(std::runtime_error) ExplorationTree::fetchMol(const std::string& canonSMILES);
%catches(std::runtime_error) ExplorationTree::deleteSubtree(const std::string& canonSMILES);
%ignore ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret);
%include "ExplorationTree.hpp";


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

// MoplherMol wrapper
%include "std_set.i"
%template(StringSet) std::set<std::string>;
%newobject MolpherMol::copy();
%ignore MolpherMol::MolpherMol(MolpherMolecule& mol);
%ignore MolpherMol::MolpherMol(MolpherMolecule& mol, bool copy);
%ignore MolpherMol::fetchMolpherMolecule();
%ignore MolpherMol::MolpherMol();
%ignore MolpherMol::operator=(const MolpherMol&);
%include "MolpherMol.hpp"

// ExplorationParameters wrapper
%template(StringVector) std::vector<std::string>;
%newobject ExplorationParameters::getSourceMol();
%newobject ExplorationParameters::getTargetMol();
%newobject ExplorationParameters::getChemOperators();
%catches(std::runtime_error) ExplorationParameters::getTargetMol();
%catches(std::runtime_error) ExplorationParameters::getChemOperators();
%catches(std::runtime_error) ExplorationParameters::getFingerprint();
%catches(std::runtime_error) ExplorationParameters::getSimilarityCoef();
%catches(std::runtime_error) ExplorationParameters::setChemOperators(const std::vector<std::string>& choices);
%catches(std::runtime_error) ExplorationParameters::setFingerprint(const std::string& fp);
%catches(std::runtime_error) ExplorationParameters::setSimilarityCoef(const std::string& coef);
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
%newobject ExplorationTree::fetchLeaves();
%newobject ExplorationTree::getParams();
%newobject ExplorationTree::fetchMol(const std::string& canonSMILES);
%newobject ExplorationTree::getCandidateMorphs();
%newobject ExplorationTree::getCandidateMorphsMask();
%catches(std::runtime_error) ExplorationTree::setCandidateMorphsMask(const std::vector<bool>&);
%catches(std::runtime_error) ExplorationTree::fetchMol(const std::string& canonSMILES);
%catches(std::runtime_error) ExplorationTree::deleteSubtree(const std::string& canonSMILES);
%catches(std::runtime_error) ExplorationTree::setParams(ExplorationParameters& params);
%ignore ExplorationTree::fetchLeaves(std::vector<MolpherMol>& ret);
%include "ExplorationTree.hpp";

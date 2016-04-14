
// MoplherMol wrapper
%ignore MolpherMol::operator=(const MolpherMol&);
%include "MolpherMol.hpp"

// ExplorationData wrapper
%catches(std::runtime_error) ExplorationData::load(const std::string& file);
%catches(std::runtime_error) ExplorationData::save(const std::string& file);
%catches(std::runtime_error) ExplorationData::setSource(const MolpherMol& mol);
%catches(std::runtime_error) ExplorationData::setCandidatesMask(const std::vector<bool>& mask);
%include "ExplorationData.hpp"

// ExplorationTree wrapper
%catches(std::runtime_error) ExplorationTree::update(const ExplorationData& data);
%catches(std::runtime_error) ExplorationTree::deleteSubtree(const std::string& canonSMILES, bool descendents_only);
%include "ExplorationTree.hpp";

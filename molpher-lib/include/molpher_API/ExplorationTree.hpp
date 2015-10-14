/* 
 * File:   ExplorationTree.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 3:10 PM
 */

#ifndef EXPLORATIONTREE_HPP
#define	EXPLORATIONTREE_HPP

#include "core/PathFinderContext.h"
#include "ExplorationTreeSnapshot.hpp"
#include "ExplorationParameters.hpp"
#include "MolpherMol.hpp"

class ExplorationTree {
    
public:
    typedef tbb::concurrent_vector<MolpherMolecule> MoleculeVector;
    typedef std::vector<bool> BoolVector;
//    typedef tbb::concurrent_vector<std::string> SmileVector;
//    typedef tbb::concurrent_hash_map<std::string, bool /*dummy*/> SmileSet;
    
private:
    PathFinderContext context;
    MoleculeVector putativeLeaves;
    BoolVector putativeLeavesMask;
    
public:
    static ExplorationTree createFromSnapshot(ExplorationTreeSnapshot snapshot);
    ExplorationTree(const std::string& sourceMolAsSMILES);
    ExplorationTree(const std::string& sourceMolAsSMILES, const ExplorationParameters& params);
    void setParams(ExplorationParameters params);
    ExplorationTreeSnapshot createSnapshot();
    
    std::vector<MolpherMol> fetchLeaves();
    MolpherMol fetchMol(const std::string& canonSMILES);
    void putativeExtend();
    void extend();

};

#endif	/* EXPLORATIONTREE_HPP */


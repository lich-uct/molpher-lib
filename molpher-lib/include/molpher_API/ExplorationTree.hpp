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
    
    friend class TreeOperation;
    
public:
    typedef tbb::concurrent_vector<MolpherMolecule> MoleculeVector;
    typedef std::vector<bool> BoolVector;
    typedef tbb::concurrent_hash_map<std::string, bool /*dummy*/> SmileSet;
//    typedef tbb::concurrent_vector<std::string> SmileVector;
    
private:
    PathFinderContext context;
    int threadCount;
    MoleculeVector candidateMoprhs;
    BoolVector candidateMorphsMask;
    ExplorationTree(IterationSnapshot& snp);
    void treeInit(IterationSnapshot& snp);
    void fetchLeaves(ExplorationTree::MoleculeVector&);
    
public:
    static ExplorationTree createFromSnapshot(ExplorationTreeSnapshot& snapshot);
    
    ExplorationTree(const std::string& sourceMolAsSMILES);
    ExplorationTree(ExplorationParameters& params);
    
    void setParams(ExplorationParameters& params);
    ExplorationTreeSnapshot createSnapshot() const;
    
    void fetchLeaves(std::vector<MolpherMol>&);
    std::vector<MolpherMol> fetchLeaves();
    void generateMorphs();
//    void extend();
    MolpherMol fetchMol(const std::string& canonSMILES);
    
    void setThreadCount(int threadCnt);
    int getThreadCount();
    std::vector<MolpherMol> getCandidateMorphs();

};

#endif	/* EXPLORATIONTREE_HPP */


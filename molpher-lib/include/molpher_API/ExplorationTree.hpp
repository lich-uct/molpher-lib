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

class TreeOperation; // forward declaration to resolve circular dependency

class ExplorationTree {
    
    friend class TreeOperation;
    
public:
    typedef tbb::concurrent_vector<MolpherMolecule> MoleculeVector;
    typedef tbb::concurrent_vector<MolpherMolecule*> MoleculePointerVector;
    typedef std::vector<bool> BoolVector;
    typedef tbb::concurrent_hash_map<std::string, bool /*dummy*/> SmileSet;
    typedef tbb::concurrent_vector<std::string> SmileVector;
    
private:
    PathFinderContext context;
    int threadCount;
    MoleculeVector candidateMoprhs;
    BoolVector candidateMorphsMask;
    ExplorationTree(IterationSnapshot& snp);
    void treeInit(IterationSnapshot& snp);
    void fetchLeaves(ExplorationTree::MoleculePointerVector&);
    
public:
    static ExplorationTree* createFromSnapshot(ExplorationTreeSnapshot& snapshot);
    
    ExplorationTree(const std::string& sourceMolAsSMILES);
    ExplorationTree(ExplorationParameters& params);
    
    void setParams(ExplorationParameters& params);
    ExplorationTreeSnapshot* createSnapshot() const;
    
    void runOperation(TreeOperation& operation);
    
    void fetchLeaves(std::vector<MolpherMol>&);
    std::vector<MolpherMol> fetchLeaves();
    MolpherMol* fetchMol(const std::string& canonSMILES);
    bool hasMol(const std::string& canonSMILES);
    void generateMorphs();
    void sortMorphs();
    void filterMorphs();
    void filterMorphs(int filters);
    void extend();
    void prune();
    
    void setThreadCount(int threadCnt);
    int getThreadCount();
    std::vector<MolpherMol> getCandidateMorphs();
    std::vector<bool> getCandidateMorphsMask(); // TODO add a bitset version, return a pointer 
    void setCandidateMorphsMask(std::vector<bool>); // TODO add a bitset version, take a reference

};

#endif	/* EXPLORATIONTREE_HPP */


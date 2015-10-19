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
    MoleculeVector putativeLeaves;
    BoolVector putativeLeavesMask;
    ExplorationTree(IterationSnapshot& snp);
    void initCandidates(IterationSnapshot& snp);
    
public:
    static ExplorationTree createFromSnapshot(ExplorationTreeSnapshot& snapshot);
    ExplorationTree(const std::string& sourceMolAsSMILES);
    ExplorationTree(ExplorationParameters& params);
    void setParams(ExplorationParameters& params);
    ExplorationTreeSnapshot createSnapshot() const;
    
    std::vector<MolpherMol> fetchLeaves();
    ExplorationTree::MoleculeVector fetchLeavesAsTBBVector();
    void putativeExtend();
//    void extend();
//    MolpherMol fetchMol(const std::string& canonSMILES);
    
    void setThreadCount(int threadCnt);
    int getThreadCount();

};

#endif	/* EXPLORATIONTREE_HPP */


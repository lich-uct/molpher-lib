/* 
 * File:   ExplorationTree.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 3:10 PM
 */

#ifndef EXPLORATIONTREE_HPP
#define	EXPLORATIONTREE_HPP

#include <memory>

#include "ExplorationData.hpp"

class TreeOperation; // forward declaration to resolve circular dependency
class FindLeavesOper;

class ExplorationTree : public std::enable_shared_from_this<ExplorationTree> {
    
    friend class TreeOperation;
    friend class FindLeavesOper;
    
public:
//    typedef tbb::concurrent_vector<MolpherMolecule> MoleculeVector;
//    typedef tbb::concurrent_vector<MolpherMolecule*> MoleculePointerVector;
//    typedef std::vector<bool> BoolVector;
//    typedef tbb::concurrent_hash_map<std::string, bool /*dummy*/> SmileSet;
//    typedef tbb::concurrent_vector<std::string> SmileVector;
    
    class ExplorationTreeImpl;
    
private:
//    PathFinderContext context;
//    int threadCount;
//    MoleculeVector candidateMoprhs;
//    BoolVector candidateMorphsMask;
//    ExplorationTree(IterationSnapshot& snp);
//    void treeInit(IterationSnapshot& snp);
//    
    std::shared_ptr<ExplorationTreeImpl> pimpl;
    
private:
    ExplorationTree();
//    ExplorationTree(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES);
    
public:
    
//    ExplorationTree(const std::string& sourceMolAsSMILES);
    //    ExplorationTree(ExplorationParameters& params);
    
    static std::shared_ptr<ExplorationTree> create(const ExplorationData& data);
//    static std::shared_ptr<ExplorationTree> create(const std::string& sourceMolAsSMILES);
    static std::shared_ptr<ExplorationTree> create(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES);
    
    std::shared_ptr<ExplorationData> asData() const;
    void update(const ExplorationData& data);
    
    void runOperation(TreeOperation& operation);
    
    std::vector<std::shared_ptr<MolpherMol> > fetchLeaves(bool increase_dist_improve_counter = false);
    std::shared_ptr<MolpherMol> fetchMol(const std::string& canonSMILES);
    bool hasMol(const std::string& canonSMILES);
    bool hasMol(std::shared_ptr<MolpherMol> mol);
//    bool isPathFound();
//    void deleteSubtree(const std::string& canonSMILES);
//    void generateMorphs();
//    void sortMorphs();
//    void filterMorphs();
//    void filterMorphs(int filters);
//    void extend();
//    void prune();
//    
//    int getThreadCount();
//    int getGenerationCount();
//    ExplorationParameters& getParams();
//    const std::vector<MolpherMol>& getCandidateMorphs();
//    const std::vector<bool>& getCandidateMorphsMask(); // TODO add a bitset version 
//    
//    void setThreadCount(int threadCnt);
//    void setParams(ExplorationParameters& params);
//    void setCandidateMorphsMask(const std::vector<bool>&); // TODO add a bitset version

};

#endif	/* EXPLORATIONTREE_HPP */


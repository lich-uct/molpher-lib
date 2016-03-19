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
#include "operations/callbacks/TraverseCallback.hpp"
#include "operations/FilterMorphsOper.hpp"

class ExplorationTree 
#ifndef SWIG
: public std::enable_shared_from_this<ExplorationTree>
#endif
 {
    
    friend class TreeOperation;
    friend class FindLeavesOper;
    friend class GenerateMorphsOper;
    friend class SortMorphsOper;
    friend class FilterMorphsOper;
    friend class ExtendTreeOper;
    friend class PruneTreeOper;
    friend class TraverseOper;
    
public:
    class ExplorationTreeImpl;
    
private:
    std::shared_ptr<ExplorationTreeImpl> pimpl;
    ExplorationTree();
    
public:
    
    static std::shared_ptr<ExplorationTree> create(const ExplorationData& data);
    static std::shared_ptr<ExplorationTree> create(const std::string& filename);
    static std::shared_ptr<ExplorationTree> create(const std::string& sourceMolAsSMILES, const std::string& targetMolAsSMILES);
    
    std::shared_ptr<ExplorationData> asData() const;
    void update(const ExplorationData& data);
    
    void runOperation(TreeOperation& operation);
    
    std::vector<std::shared_ptr<MolpherMol> > fetchLeaves(bool increase_dist_improve_counter = false);
    std::shared_ptr<MolpherMol> fetchMol(const std::string& canonSMILES);
    bool hasMol(const std::string& canonSMILES);
    bool hasMol(std::shared_ptr<MolpherMol> mol);
    bool isPathFound();
    void deleteSubtree(const std::string& canonSMILES, bool descendents_only = false);
    void generateMorphs();
    void sortMorphs();
    void filterMorphs(bool verbose_output = false);
    void filterMorphs(FilterMorphsOper::MorphFilters filters, bool verbose_output = false);
    void extend();
    void prune();
    void traverse(const std::string& rootSMILES, TraverseCallback& callback);
    void traverse(TraverseCallback& callback);
    void save(const std::string& filename);
//    
//    int getThreadCount();
    unsigned getGenerationCount();
    std::vector<std::shared_ptr<MolpherMol> > getCandidateMorphs();
    std::vector<bool> getCandidateMorphsMask(); // TODO add a bitset version 
//    
//    void setThreadCount(int threadCnt);
//    void setCandidateMorphsMask(const std::vector<bool>&); // TODO add a bitset version

};

#endif	/* EXPLORATIONTREE_HPP */


/* 
 * File:   ExplorationData.hpp
 * Author: sichom
 *
 * Created on February 9, 2016, 12:40 PM
 */

#ifndef EXPLORATIONDATA_HPP
#define	EXPLORATIONDATA_HPP

#include <string>
#include <vector>
#include <map>

#include "data_structs/MolpherMol.hpp"

class ExplorationData {    
    
private:
    struct ExplorationDataImpl;
    std::unique_ptr<ExplorationDataImpl> pimpl;

public:
    
    ExplorationData();
    
    // getters
    
    unsigned getGenerationCount() const;
    unsigned getThreadCount() const;
    int getFingerprint() const;
    int getSimilarityCoefficient() const;
    std::set<int> getChemicalOperators() const;
    
    // params
    double getMinAcceptableMolecularWeight() const; // weightMin
    double getMaxAcceptableMolecularWeight() const; // weightMax
    int getCntCandidatesToKeep() const; // acceptMin
    int getCntCandidatesToKeepMax() const; // acceptMax
    int getCntMorphs() const; // farProduce
    int getCntMorphsInDepth() const; // closeProduce
    double getDistToTargetDepthSwitch() const; // farCloseThreshold (if below this threshold, generate closeProduce morphs)
    int getCntMaxMorphs() const; // maxMorhpsTotal (prune this and all descendents if no distance improvement among children for this many iters)
    int getItThreshold() const; // nonProducingSurvive (prune descendents if no distance improvement among children after this many iters)
    
    std::unique_ptr<MolpherMol> getSource() const;
    std::unique_ptr<MolpherMol> getTarget() const;

    std::unique_ptr<std::vector<std::unique_ptr<MolpherMol> > > getCandidates() const;
    const std::vector<bool>& getCandidatesMask() const;
    std::unique_ptr<std::map<std::string, std::unique_ptr<MolpherMol> > > getTreeMap() const;
    const std::map<std::string, unsigned>& getDerivationMap() const;
    
    // setters
    
    void setGenerationCount(unsigned);
    void setThreadCount(unsigned);
    void setFingerprint(int);
    void setSimilarityCoefficient(int);
    void setChemicalOperators(const std::set<int>&);
    void addChemicalOperator(int);
    void removeChemicalOperator(int);
    
    // params
    void setMinAcceptableMolecularWeight(double); // weightMin
    void setMaxAcceptableMolecularWeight(double); // weightMax
    void setCntCandidatesToKeep(int); // acceptMin
    void setCntCandidatesToKeepMax(int); // acceptMax
    void setCntMorphs(int); // farProduce
    void setCntMorphsInDepth(int); // closeProduce
    void setDistToTargetDepthSwitch(double); // farCloseThreshold (if below this threshold, generate closeProduce morphs)
    void setCntMaxMorphs(int); // maxMorhpsTotal (prune this and all descendents if no distance improvement among children for this many iters)
    void setItThreshold(int); // nonProducingSurvive (prune descendents if no distance improvement among children after this many iters)
    
    void setSource(const MolpherMol&);
    void setTarget(const MolpherMol&);

    void setCandidates(const std::vector<MolpherMol>&);
    void addCandidate(const MolpherMol&, unsigned index);
    void addCandidate(const MolpherMol&);
    void removeCandidate(unsigned index);
    void setCandidatesMask(const std::vector<bool>&);
    void setCandidatesMaskAt(bool, unsigned);
//    void setTreeMap(const std::map<std::string, MolpherMol>&);
    void addToTreeMap(const std::string&, const MolpherMol&);
    std::unique_ptr<MolpherMol> popFromTreeMap(const std::string&);
//    void setDerivationMap(const std::map<std::string, unsigned>&);
    void addToDerivationMap(const std::string&, unsigned);
    void increaseDerivationsCount(const std::string&);
    void decreaseDerivationsCount(const std::string&);
    unsigned popFromDerivationMap(const std::string&);
    
    bool isValid() const;
};

#endif	/* EXPLORATIONDATA_HPP */


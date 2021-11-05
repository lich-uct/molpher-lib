/*
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EXPLORATIONDATA_HPP
#define	EXPLORATIONDATA_HPP

#include <string>
#include <vector>
#include <map>

#include "data_structs/MolpherMol.hpp"

class ExplorationData {

    // TODO: add advanced integrity checking (if this instance is valid,
    // it should always represent a consistent tree)
    
public:
    struct ExplorationDataImpl;
    
private:
    std::unique_ptr<ExplorationDataImpl> pimpl;

public:
    
    ExplorationData();
    ~ExplorationData();
    
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
	double getSAScoreMax() const;
    
    std::shared_ptr<MolpherMol> getSource() const;
    std::shared_ptr<MolpherMol> getTarget() const;

    std::shared_ptr<std::vector<std::shared_ptr<MolpherMol> > > getCandidates() const;
    const std::vector<bool>& getCandidatesMask() const;
    std::shared_ptr<std::map<std::string, std::shared_ptr<MolpherMol> > > getTreeMap() const;
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
	void setSAScoreMax(double);

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
    std::shared_ptr<MolpherMol> popFromTreeMap(const std::string&);
//    void setDerivationMap(const std::map<std::string, unsigned>&);
    void addToDerivationMap(const std::string&, unsigned);
    void increaseDerivationsCount(const std::string&);
    void decreaseDerivationsCount(const std::string&);
    unsigned popFromDerivationMap(const std::string&);
    
    static std::shared_ptr<ExplorationData> load(const std::string &file);
    void save(const std::string &file);
    
    bool isValid() const;
};

#endif	/* EXPLORATIONDATA_HPP */


/* 
 * File:   ExplorationParameters.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 12:31 PM
 */

#ifndef EXPLORATIONPARAMETERS_HPP
#define	EXPLORATIONPARAMETERS_HPP

#include "IterationSnapshot.h"
#include "MolpherMol.hpp"
#include "chemoper_selectors.h"
#include "simcoeff_selectors.h"
#include "fingerprint_selectors.h"

class ExplorationParameters {
    
    friend class ExplorationTree;
    
private:
    IterationSnapshot iterSnapshot;
    static std::map<std::string, ChemOperSelector> CHEMOPERS_SELECTOR_MAPPPING;
    static std::map<std::string, SimCoeffSelector> SIMILARITY_SELECTOR_MAPPPING;
    static std::map<std::string, FingerprintSelector> FINGERPRINT_SELECTOR_MAPPPING;
    
public:
    ExplorationParameters();
//    IterationSnapshot createIterationSnapshot() const;
    bool valid();
    
    MolpherMol* getSourceMol();
    MolpherMol* getTargetMol();
    const std::vector<std::string>& getChemOperators();
    std::string getFingerprint();
    std::string getSimilarityCoef();
    double getMinAcceptableMolecularWeight(); // weightMin
    double getMaxAcceptableMolecularWeight(); // weightMax
    int getCntCandidatesToKeep(); // acceptMin
    int getCntCandidatesToKeepMax(); // acceptMax
    int getCntMorphs(); // farProduce
    int getCntMorphsInDepth(); // closeProduce
    double getDistToTargetDepthSwitch(); // farCloseThreashold (if below this threshold, generate closeProduce morphs)
    int getCntMaxMorphs(); // maxMorhpsTotal (prune this and all descendents if no distance improvement among children for this many iters)
    int getItThreshold(); // nonProducingSurvive (prune descendents if no distance improvement among children after this many iters)
    
    void setSourceMol(const std::string& mol);
    void setSourceMol(MolpherMol& mol);
    void setTargetMol(const std::string& mol);
    void setTargetMol(MolpherMol& mol);
    void setChemOperators(const std::vector<std::string>& choices);
    void setFingerprint(const std::string& fp);
    void setSimilarityCoef(const std::string& coef);
    void setMinAcceptableMolecularWeight(double weight); // weightMin
    void setMaxAcceptableMolecularWeight(double weight); // weightMax
    void setCntCandidatesToKeep(int value); // acceptMin
    void setCntCandidatesToKeepMax(int value); // acceptMax
    void setCntMorphs(int value); // farProduce
    void setCntMorphsInDepth(int value); // closeProduce
    void setDistToTargetDepthSwitch(double value); // farCloseThreashold (if below this threshold, generate closeProduce morphs)
    void setCntMaxMorphs(int value); // maxMorhpsTotal (prune this and all descendents if no distance improvement among children for this many iters)
    void setItThreshold(int value); // nonProducingSurvive (prune descendents if no distance improvement among children after this many iters)
    
};

#endif	/* EXPLORATIONPARAMETERS_HPP */


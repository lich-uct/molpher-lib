/* 
 * File:   MolpherMol.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 1:25 PM
 */

#ifndef MOLPHERMOL_HPP
#define	MOLPHERMOL_HPP

#include <string>
#include <memory>
#include <set>

//#include "MolpherMolecule.h"

#include "selectors/chemoper_selectors.h"

class ExplorationTree; // forward declaration to resolve circular dependency

class MolpherMol {
    
private:
//    MolpherMolecule* mol;
//    bool selfAllocated;
    class MolpherMolImpl;
    std::unique_ptr<MolpherMolImpl> pimpl;

public:
//    MolpherMol(std::shared_ptr<MolpherMolImpl> pimpl);
    
    MolpherMol();
//    MolpherMol(MolpherMolecule& mol);
//    MolpherMol(MolpherMolecule& mol, bool copy);
    MolpherMol(const std::string& smiles, const std::string& formula, const std::string& parentSmile,
                const unsigned& oper, const double& dist, const double& distToClosestDecoy,
                const double& weight, const double& sascore);
    MolpherMol(const std::string& smiles);
    MolpherMol(const MolpherMol& other);
    ~MolpherMol();
    
    MolpherMol& operator=(const MolpherMol&);
    
//    MolpherMolecule& fetchMolpherMolecule() const; // TODO: maybe get rid of this
    std::shared_ptr<MolpherMol> copy() const;
    
    // getters
    const std::string& getSMILES() const;
    double getDistToTarget() const;
    std::shared_ptr<ExplorationTree> getTree();
    const std::string& getParentSMILES() const;
    const std::set<std::string>& getDescendants() const;
    const std::set<std::string>& getHistoricDescendants() const;
    unsigned int getItersWithoutDistImprovement() const;
    double getSAScore() const;
    double getMolecularWeight() const;
    const std::string& getFormula() const;
    int getParentOper() const;
    
    // setters
    void setOwner(std::shared_ptr<ExplorationTree> tree);
    void setSMILES(const std::string&);
    void setParentSMILES(const std::string&);
    void setDistToTarget(double dist);
    void setSAScore(double score);
    void setItersWithoutDistImprovement(unsigned int count);
    void increaseItersWithoutDistImprovement();
    void decreaseItersWithoutDistImprovement();
    void addToDescendants(const std::string& smiles);
    void removeFromDescendants(const std::string& smiles);
    void setDescendants(const std::set<std::string>&);
    void addToHistoricDescendants(const std::string& smiles);
    void removeFromHistoricDescendants(const std::string& smiles);
    void setHistoricDescendants(const std::set<std::string>&);
    
    bool isValid() const;
    bool isBoundToTree() const;
};

#endif	/* MOLPHERMOL_HPP */


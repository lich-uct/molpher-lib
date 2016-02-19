/* 
 * File:   MolpherMolImpl.hpp
 * Author: sichom
 *
 * Created on January 27, 2016, 3:47 PM
 */

#ifndef MOLPHERMOLIMPL_HPP
#define	MOLPHERMOLIMPL_HPP

#include <memory>
#include <vector>

#include "core/data_structs/MolpherMolData.hpp"
#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"

class MolpherMol::MolpherMolImpl {
    
public:
    std::shared_ptr<ExplorationTree> tree;
    MolpherMolData data;
    
    MolpherMolImpl(const std::string& SMILES);
    MolpherMolImpl(const MolpherMolData& data);
    MolpherMolImpl(const MolpherMolImpl& other);
    MolpherMolImpl();
    
    std::unique_ptr<MolpherMolImpl> copy() const;
//    MolpherMolData asData() const;
//     
//    std::string getSMILES() const;
//    double getDistToTarget() const;
//    std::shared_ptr<ExplorationTree> getTree();
//    std::string getParentSMILES() const;
//    const std::unique_ptr<std::set<std::string> >& getDescendants() const;
//    const std::unique_ptr<std::set<std::string> >& getHistoricDescendants() const;
//    unsigned int getItersWithoutDistImprovement() const;
//    double getSAScore() const;
//    double getMolecularWeight() const;
//    
    void setSMILES(const std::string& smiles);
//    void setDistToTarget(double dist);
//    void setSAScore(double score);
//    void setItersWithoutDistImprovement(unsigned int count);
//    void increaseItersWithoutDistImprovement();
//    void decreaseItersWithoutDistImprovement();
//    void addToDescendants(const std::string& smiles);
//    void removeFromDescendants(const std::string& smiles);
//    void addToHistoricDescendants(const std::string& smiles);
//    void removeFromHistoricDescendants(const std::string& smiles);
//    
//    bool isValid() const;
//    bool isBoundToTree() const;
};

#endif	/* MOLPHERMOLIMPL_HPP */


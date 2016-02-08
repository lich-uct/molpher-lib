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
    
private:
    std::shared_ptr<ExplorationTree> tree;
    MolpherMolData data;
    
public:
    MolpherMolImpl(const std::string& SMILES);
    MolpherMolImpl(const MolpherMolData& data);
    MolpherMolImpl(const MolpherMolImpl& other);
    MolpherMolImpl();
    
    std::unique_ptr<MolpherMolImpl> copy() const;
    MolpherMolData asData() const;
//    bool isValid() const;
//    
    std::string getSMILES() const;
    double getDistToTarget() const;
    std::shared_ptr<ExplorationTree> getTree();
//    std::string getParentSMILES() const;
//    std::shared_ptr<std::vector<std::shared_ptr<MolpherMolImpl> > > getDescendants() const;
//    std::shared_ptr<std::vector<std::shared_ptr<MolpherMolImpl> > > getHistoricDescendants() const;
//    unsigned int getItersWithoutDistImprovement() const;
//    double getSAScore() const;
//    double getMolecularWeight() const;
//    
    void setDistToTarget(double dist);
//    void setSAScore(double score);
//    void setItersWithoutDistImprovement(unsigned int count);
    
};

#endif	/* MOLPHERMOLIMPL_HPP */


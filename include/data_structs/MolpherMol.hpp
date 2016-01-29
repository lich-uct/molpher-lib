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

//#include "MolpherMolecule.h"

class MolpherMol {
    
private:
//    MolpherMolecule* mol;
//    bool selfAllocated;
    
    std::shared_ptr<MolpherMolImpl> pimpl;
    
public:
    class MolpherMolImpl;
    
    MolpherMol(std::shared_ptr<MolpherMolImpl> pimpl);
//    MolpherMol();
//    MolpherMol(MolpherMolecule& mol);
//    MolpherMol(MolpherMolecule& mol, bool copy);
//    MolpherMol(const MolpherMol& other);
//    ~MolpherMol();
    
//    MolpherMol& operator=(const MolpherMol&);
//    
//    MolpherMolecule& fetchMolpherMolecule() const; // TODO: maybe get rid of this
//    bool isBound() const;
//    MolpherMol* copy() const;
//    
//    std::string getSMILES();
//    double getDistToTarget();
//    std::string getParentSMILES();
//    const std::set<std::string>& getDescendants();
//    const std::set<std::string>& getHistoricDescendants();
//    unsigned int getItersWithoutDistImprovement();
//    double getSAScore();
//    double getMolecularWeight();
//    
//    void setDistToTarget(double dist);
//    void setSAScore(double score);
};

#endif	/* MOLPHERMOL_HPP */


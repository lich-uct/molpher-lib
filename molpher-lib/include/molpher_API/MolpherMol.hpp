/* 
 * File:   MolpherMol.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 1:25 PM
 */

#ifndef MOLPHERMOL_HPP
#define	MOLPHERMOL_HPP

#include "MolpherMolecule.h"
#include <string>

class MolpherMol {
    
private:
    MolpherMolecule* mol;
    bool selfAllocated;
    
public:
    MolpherMol();
    MolpherMol(MolpherMolecule& mol);
    MolpherMol(MolpherMolecule& mol, bool copy);
    MolpherMol(const MolpherMol& other);
    ~MolpherMol();
    
    MolpherMol& operator=(const MolpherMol&);
    
    MolpherMolecule& fetchMolpherMolecule() const; // TODO: maybe get rid of this
    bool isBound() const;
    MolpherMol* copy() const;
    
    std::string getSMILES();
    double getDistToTarget();
    std::string getParentSMILES();
    const std::set<std::string>& getDescendants();
    const std::set<std::string>& getHistoricDescendants();
    unsigned int getItersWithoutDistImprovement();
    double getSAScore();
    double getMolecularWeight();
    
    void setDistToTarget(double dist);
    void setSAScore(double dist);
};

#endif	/* MOLPHERMOL_HPP */


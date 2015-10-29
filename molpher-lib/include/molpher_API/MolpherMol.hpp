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
    // TODO: memory management
    
private:
    MolpherMolecule* mol;
    
public:
    MolpherMol();
    MolpherMol(MolpherMolecule& mol);
    MolpherMol(const std::string &smile);
    MolpherMolecule& getMol(); // TODO: get rid of this
    
    std::string getSMILES();
    double getDistToTarget();
    std::string getParentSMILES();
    // TODO create more getters for the members of mol
};

#endif	/* MOLPHERMOL_HPP */


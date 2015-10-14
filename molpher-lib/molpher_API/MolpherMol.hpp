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
    MolpherMolecule mol;
    
public:
    MolpherMol(const std::string &smile);
    MolpherMolecule& getMol();
    
    // TODO create getters for the members of mol
};

#endif	/* MOLPHERMOL_HPP */


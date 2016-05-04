/* 
 * File:   MolpherAtom.h
 * Author: Petyr
 *
 * Created on 3.04.2013, 12:34
 */

#pragma once

#include "GraphMol/Atom.h"

// used type definition
typedef int AtomicNum;

/**
 * Molpher atom representation.
 */
struct MolpherAtom
{
public:    
    MolpherAtom(RDKit::Atom *atom) :
        atomicNum(atom->getAtomicNum()), formalCharge(atom->getFormalCharge()),
        mass(atom->getMass())
    { }
public:
    
    AtomicNum atomicNum;
    
    int formalCharge;
    
    double mass;
};



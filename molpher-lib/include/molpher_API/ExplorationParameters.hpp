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

class ExplorationParameters {
    
    friend class ExplorationTree;
    
private:
    IterationSnapshot iterSnapshot;
    
public:
    ExplorationParameters();
//    IterationSnapshot createIterationSnapshot() const;
    bool valid();
    
    void setSourceMol(const std::string& mol);
    void setSourceMol(MolpherMol& mol);
    MolpherMol* getSourceMol() const;
    // TODO more getters and setters for the parameters
    
};

#endif	/* EXPLORATIONPARAMETERS_HPP */


/* 
 * File:   ExplorationParameters.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 12:31 PM
 */

#ifndef EXPLORATIONPARAMETERS_HPP
#define	EXPLORATIONPARAMETERS_HPP

#include "IterationSnapshot.h"

class ExplorationParameters {
private:
    IterationSnapshot iterSnapshot;
    
public:
    ExplorationParameters();
    IterationSnapshot createIterationSnapshot();
    
    // TODO getters and setters for the parameters
    
};

#endif	/* EXPLORATIONPARAMETERS_HPP */


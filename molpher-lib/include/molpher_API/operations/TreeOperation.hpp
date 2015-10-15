/* 
 * File:   TreeOperation.hpp
 * Author: sichom
 *
 * Created on October 13, 2015, 12:43 PM
 */

#ifndef TREEOPERATION_HPP
#define	TREEOPERATION_HPP

#include "molpher_API/ExplorationTree.hpp"
#include <string>

class TreeOperation {

protected:
    ExplorationTree& tree;
    const int threadCnt;
    PathFinderContext& fetchTreeContext();
    
public:
    TreeOperation(ExplorationTree& expTree);
    virtual void operator()() = 0;
};

#endif	/* TREEOPERATION_HPP */


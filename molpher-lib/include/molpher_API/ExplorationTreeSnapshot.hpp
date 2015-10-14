/* 
 * File:   ExplorationTreeSnapshot.hpp
 * Author: sichom
 *
 * Created on October 14, 2015, 2:09 PM
 */

#ifndef EXPLORATIONTREESNAPSHOT_HPP
#define	EXPLORATIONTREESNAPSHOT_HPP

#include "IterationSnapshot.h"

class ExplorationTreeSnapshot {

private:
    IterationSnapshot iterSnapshot;
    
public:
    static ExplorationTreeSnapshot load(const std::string& filename);
    void save(const std::string& filename);
    ExplorationTreeSnapshot(IterationSnapshot iterSnapshot);
};

#endif	/* EXPLORATIONTREESNAPSHOT_HPP */


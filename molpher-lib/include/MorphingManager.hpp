/* 
 * File:   MorphingManager.hpp
 * Author: sichom
 *
 * Created on October 7, 2015, 12:37 PM
 */

#ifndef MORPHINGMANAGER_HPP
#define	MORPHINGMANAGER_HPP

#include "IterationSnapshot.h"
#include "core/PathFinderContext.h"

class MoprhingManager {
private:

    IterationSnapshot mSnapshot;
    PathFinderContext mContext;
    int mThreadCount;
    
    void initContext();

public:
    //        MoprhingManager();
    MoprhingManager(const std::string &jobFile);
    MoprhingManager(const std::string &jobFile, int threadCnt);
};

#endif	/* MORPHINGMANAGER_HPP */


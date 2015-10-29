
#ifndef TRAVERSECALLBACK_HPP
#define	TRAVERSECALLBACK_HPP

#include "core/PathFinderContext.h"

#include "../MolpherMol.hpp"

class TraverseCallback {
    // TODO: add a method to the tree that takes this and runs TraverseOper with it and test calling from python
    
    protected:
        PathFinderContext* context;
    
    public:
        TraverseCallback(PathFinderContext& context);
        TraverseCallback();
        virtual ~TraverseCallback();
        virtual void processMorph(MolpherMol& morph) = 0;
};

#endif	/* TRAVERSECALLBACK_HPP */
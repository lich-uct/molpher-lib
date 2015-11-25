
#ifndef TRAVERSECALLBACK_HPP
#define	TRAVERSECALLBACK_HPP

#include "core/PathFinderContext.h"

#include "../MolpherMol.hpp"

class TraverseCallback {
    
    protected:
        PathFinderContext* context;
    
    public:
        TraverseCallback(PathFinderContext& context);
        TraverseCallback();
        virtual ~TraverseCallback();
        virtual void processMorph(MolpherMol& morph) = 0;
};

#endif	/* TRAVERSECALLBACK_HPP */
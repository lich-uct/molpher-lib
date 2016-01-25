
#ifndef TRAVERSECALLBACK_HPP
#define	TRAVERSECALLBACK_HPP

#include "core/misc/PathFinderContext.h"

#include "data_structs/MolpherMol.hpp"

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
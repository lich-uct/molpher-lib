
#ifndef TRAVERSECALLBACK_HPP
#define	TRAVERSECALLBACK_HPP

#include "data_structs/MolpherMol.hpp"

class TraverseCallback {    
    public:
        class TraverseCallbackImpl;
        
        TraverseCallback();
        virtual ~TraverseCallback();
        virtual void operator()(std::shared_ptr<MolpherMol> morph) const = 0;
        
    protected:
        void setTraverseCallbackPimpl(std::shared_ptr<TraverseCallbackImpl> pimpl);
        
    private:
        std::shared_ptr<TraverseCallbackImpl> pimpl;
};

#endif	/* TRAVERSECALLBACK_HPP */
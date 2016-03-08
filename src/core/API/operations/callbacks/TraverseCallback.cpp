
#include "operations/callbacks/TraverseCallback.hpp"
#include "TraverseCallbackImpl.hpp"

TraverseCallback::TraverseCallback() : 
pimpl(new TraverseCallback::TraverseCallbackImpl())
{
    // no action
}


TraverseCallback::~TraverseCallback() {
    // no action
}

void TraverseCallback::setTraverseCallbackPimpl(std::shared_ptr<TraverseCallbackImpl> pimpl) {
    this->pimpl = pimpl;
}

// pimpl

TraverseCallback::TraverseCallbackImpl::TraverseCallbackImpl() {
    // no action
}

//void TraverseCallback::TraverseCallbackImpl::operator()(std::shared_ptr<MolpherMol> morph) {
//    throw std::runtime_error("This is just a dummy method. It should not be called directly!");
//}


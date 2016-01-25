
%{
#include "operations/callbacks/TraverseCallback.hpp"
%}

%ignore TraverseCallback::TraverseCallback(PathFinderContext& context);
%include "TraverseCallback.hpp"
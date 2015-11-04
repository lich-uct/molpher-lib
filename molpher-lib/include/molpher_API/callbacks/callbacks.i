
%{
#include "molpher_API/callbacks/TraverseCallback.hpp"
%}

%ignore TraverseCallback::TraverseCallback(PathFinderContext& context);
%include "TraverseCallback.hpp"
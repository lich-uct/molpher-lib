
%{
#include "molpher_API/operations/TreeOperation.hpp"
#include "molpher_API/operations/FindLeavesOper.hpp"
#include "molpher_API/operations/GenerateMorphsOper.hpp"
#include "molpher_API/operations/SortMorphsOper.hpp"
#include "molpher_API/operations/FilterMorphsOper.hpp"
#include "molpher_API/operations/ExtendTreeOper.hpp"
#include "molpher_API/operations/PruneTreeOper.hpp"
#include "molpher_API/operations/TraverseOper.hpp"
#include "molpher_API/callbacks/TraverseCallback.hpp"
%}

%include "TreeOperation.hpp"
%newobject FindLeavesOper::fetchLeaves();        
%include "FindLeavesOper.hpp"
%include "GenerateMorphsOper.hpp"
%include "SortMorphsOper.hpp"
%include "FilterMorphsOper.hpp"
%include "ExtendTreeOper.hpp"
%include "PruneTreeOper.hpp"
%include "TraverseOper.hpp"

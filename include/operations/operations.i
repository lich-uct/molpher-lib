
%{
#include "operations/TreeOperation.hpp"
#include "operations/FindLeavesOper.hpp"
//#include "operations/GenerateMorphsOper.hpp"
//#include "operations/SortMorphsOper.hpp"
//#include "operations/FilterMorphsOper.hpp"
//#include "operations/ExtendTreeOper.hpp"
//#include "operations/PruneTreeOper.hpp"
//#include "operations/TraverseOper.hpp"
//#include "operations/callbacks/TraverseCallback.hpp"
%}

%template(MolVector) std::vector<std::shared_ptr<MolpherMol> >;

%include "TreeOperation.hpp"
//%newobject FindLeavesOper::fetchLeaves();
%catches(std::runtime_error) FindLeavesOper::operator()();
%include "FindLeavesOper.hpp"
//%catches(std::runtime_error) GenerateMorphsOper::operator()();
//%include "GenerateMorphsOper.hpp"
//%catches(std::runtime_error) SortMorphsOper::operator()();
//%include "SortMorphsOper.hpp"
//%catches(std::runtime_error) FilterMorphsOper::operator()();
//%include "FilterMorphsOper.hpp"
//%catches(std::runtime_error) ExtendTreeOper::operator()();
//%include "ExtendTreeOper.hpp"
//%catches(std::runtime_error) PruneTreeOper::operator()();
//%include "PruneTreeOper.hpp"
//%catches(std::runtime_error) TraverseOper::operator()();
//%include "TraverseOper.hpp"

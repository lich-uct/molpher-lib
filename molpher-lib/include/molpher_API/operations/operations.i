
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
%catches(std::runtime_error) FindLeavesOper::operator()();
%include "FindLeavesOper.hpp"
%catches(std::runtime_error) GenerateMorphsOper::operator()();
%include "GenerateMorphsOper.hpp"
%catches(std::runtime_error) SortMorphsOper::operator()();
%include "SortMorphsOper.hpp"
%catches(std::runtime_error) FilterMorphsOper::operator()();
%include "FilterMorphsOper.hpp"
%catches(std::runtime_error) ExtendTreeOper::operator()();
%include "ExtendTreeOper.hpp"
%catches(std::runtime_error) PruneTreeOper::operator()();
%include "PruneTreeOper.hpp"
%catches(std::runtime_error) TraverseOper::operator()();
%include "TraverseOper.hpp"


%ignore TreeOperation::TreeOperationImpl;
%include "TreeOperation.hpp";

%ignore FindLeavesOper::FindLeavesOperImpl;
%catches(std::runtime_error) FindLeavesOper::operator()();
%include "FindLeavesOper.hpp";

%ignore GenerateMorphsOper::GenerateMorphsOperImpl;
%catches(std::runtime_error) GenerateMorphsOper::operator()();
%include "GenerateMorphsOper.hpp";

%ignore SortMorphsOper::SortMorphsOperImpl;
%catches(std::runtime_error) SortMorphsOper::operator()();
%include "SortMorphsOper.hpp";

%ignore FilterMorphsOper::FilterMorphsOperImpl;
%catches(std::runtime_error) FilterMorphsOper::operator()();
%include "FilterMorphsOper.hpp";

%ignore ExtendTreeOper::ExtendTreeOperImpl;
%catches(std::runtime_error) ExtendTreeOper::operator()();
%include "ExtendTreeOper.hpp";

%ignore PruneTreeOper::PruneTreeOperImpl;
%catches(std::runtime_error) PruneTreeOper::operator()();
%include "PruneTreeOper.hpp";

%ignore TraverseOper::TraverseOperImpl;
%catches(std::runtime_error) TraverseOper::operator()();
%include "TraverseOper.hpp";

/*
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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

%ignore CleanMorphsOper::CleanMorphsOperImpl;
%catches(std::runtime_error) CleanMorphsOper::operator()();
%include "CleanMorphsOper.hpp";

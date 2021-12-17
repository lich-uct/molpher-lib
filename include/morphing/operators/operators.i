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

%ignore MorphingOperator::MorphingOperatorImpl;
%include "MorphingOperator.hpp";

%ignore AddAtom::AddAtomImpl;
%catches(std::runtime_error) AddAtom::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) AddAtom::morph();
%include "AddAtom.hpp";

%ignore RemoveAtom::RemoveAtomImpl;
%catches(std::runtime_error) RemoveAtom::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) RemoveAtom::morph();
%include "RemoveAtom.hpp";

%ignore MutateAtom::MutateAtomImpl;
%catches(std::runtime_error) MutateAtom::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) MutateAtom::morph();
%include "MutateAtom.hpp";

%ignore AddBond::AddBondImpl;
%catches(std::runtime_error) AddBond::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) AddBond::morph();
%include "AddBond.hpp";

%ignore RemoveBond::RemoveBondImpl;
%catches(std::runtime_error) RemoveBond::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) RemoveBond::morph();
%include "RemoveBond.hpp";

%ignore AddAtom::AddAtomImpl;
%catches(std::runtime_error) ContractBond::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) ContractBond::morph();
%include "ContractBond.hpp";

%ignore InterlayAtom::InterlayAtomImpl;
%catches(std::runtime_error) InterlayAtom::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) InterlayAtom::morph();
%include "InterlayAtom.hpp";

%ignore RerouteBond::RerouteBondImpl;
%catches(std::runtime_error) RerouteBond::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) RerouteBond::morph();
%include "RerouteBond.hpp";

%ignore ReactionOperator::ReactionOperatorImpl;
%catches(std::runtime_error) ReactionOperator::setOriginal(std::shared_ptr<MolpherMol> mol);
%catches(std::runtime_error) ReactionOperator::morph();
%include "ReactionOperator.hpp";
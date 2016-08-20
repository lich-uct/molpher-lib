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

#ifndef MOLPHER_LIB_CLEANMORPHSOPER_HPP
#define MOLPHER_LIB_CLEANMORPHSOPER_HPP

#include "TreeOperation.hpp"

class CleanMorphsOper : public TreeOperation {

public:
	class CleanMorphsOperImpl;

	CleanMorphsOper(std::shared_ptr<ExplorationTree> expTree);
	CleanMorphsOper();
	virtual void operator()();

private:
	std::shared_ptr<CleanMorphsOperImpl> pimpl;
};

#endif //MOLPHER_LIB_CLEANMORPHSOPER_HPP

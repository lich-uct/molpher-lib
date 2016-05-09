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

#ifndef PRUNETREEOPER_HPP
#define	PRUNETREEOPER_HPP

#include "TreeOperation.hpp"

class PruneTreeOper : public TreeOperation {
    
    public:
        class PruneTreeOperImpl;
        PruneTreeOper(std::shared_ptr<ExplorationTree> expTree);
        PruneTreeOper();
        virtual void operator()();
        
//        const std::vector<std::string>& getPruned();
        
    private:
        std::shared_ptr<PruneTreeOperImpl> pimpl;
        
};

#endif	/* PRUNETREEOPER_HPP */


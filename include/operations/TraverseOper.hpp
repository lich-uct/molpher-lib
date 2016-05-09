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

#ifndef TRAVERSEOPER_HPP
#define	TRAVERSEOPER_HPP

#include "TreeOperation.hpp"
#include "operations/callbacks/TraverseCallback.hpp"

class TraverseOper : public TreeOperation {
    
    public:
        class TraverseOperImpl;
        
        TraverseOper(std::shared_ptr<ExplorationTree> expTree, TraverseCallback& callback);
        TraverseOper(TraverseCallback& callback);
        TraverseOper(std::shared_ptr<ExplorationTree> expTree, TraverseCallback& callback, const std::string& rootSMILES);
        virtual void operator()();
        
    private:
        std::shared_ptr<TraverseOperImpl> pimpl;
};


#endif	/* TRAVERSEOPER_HPP */


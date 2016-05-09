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

#ifndef TRAVERSECALLBACK_HPP
#define	TRAVERSECALLBACK_HPP

#include "data_structs/MolpherMol.hpp"

class TraverseCallback {    
    public:
        class TraverseCallbackImpl;
        
        TraverseCallback();
        virtual ~TraverseCallback();
        virtual void operator()(std::shared_ptr<MolpherMol> morph) const = 0;
        
    protected:
        void setTraverseCallbackPimpl(std::shared_ptr<TraverseCallbackImpl> pimpl);
        
    private:
        std::shared_ptr<TraverseCallbackImpl> pimpl;
};

#endif	/* TRAVERSECALLBACK_HPP */
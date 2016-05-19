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

#ifndef FILTERMORPHSOPER_HPP
#define	FILTERMORPHSOPER_HPP

#include "TreeOperation.hpp"

class FilterMorphsOper : public TreeOperation {

public:
    class FilterMorphsOperImpl;

    enum MorphFilters : int {
        PROBABILITY = 1 << 0,
        WEIGHT = 1 << 1,
        SYNTHESIS = 1 << 2,
        MAX_DERIVATIONS = 1 << 3,
        DUPLICATES = 1 << 4,
        HISTORIC_DESCENDENTS = 1 << 5,
        ALL = PROBABILITY | WEIGHT | SYNTHESIS | MAX_DERIVATIONS | DUPLICATES | HISTORIC_DESCENDENTS
    };

    FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, bool verbose = false);
    FilterMorphsOper(bool verbose = false);
    FilterMorphsOper(std::shared_ptr<ExplorationTree> expTree, MorphFilters filters, bool verbose = false);
    FilterMorphsOper(MorphFilters filters, bool verbose = false);
    virtual void operator()();
    
private:
    std::shared_ptr<FilterMorphsOperImpl> pimpl;
};


#endif	/* FILTERMORPHSOPER_HPP */


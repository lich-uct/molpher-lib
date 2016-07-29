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

#ifndef MOLPHER_LIB_SORTMORPHSCALLBACK_H
#define MOLPHER_LIB_SORTMORPHSCALLBACK_H

#include "data_structs/MolpherMol.hpp"

class SortMorphsCallback {
public:
    class SortMorphsCallbackImpl;

    SortMorphsCallback();
    virtual ~SortMorphsCallback();
    virtual bool operator()(std::shared_ptr<MolpherMol> morph_1, std::shared_ptr<MolpherMol> morph_2) const = 0;

protected:
    void setSortMorphsPimpl(std::shared_ptr<SortMorphsCallbackImpl> pimpl);

private:
    std::shared_ptr<SortMorphsCallbackImpl> pimpl;
};

class DefaultSortCallback : public SortMorphsCallback {
public:
    DefaultSortCallback();
    bool operator()(std::shared_ptr<MolpherMol> morph_1, std::shared_ptr<MolpherMol> morph_2) const;
};

#endif //MOLPHER_LIB_SORTMORPHSCALLBACK_H

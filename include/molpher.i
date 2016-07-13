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

%module(directors="1", threads="1") core // FIXME: not too cool to be releasing GIL on every call (http://swig.10945.n7.nabble.com/How-to-release-Python-GIL-td5027.html)
%feature("director") TreeOperation;
%feature("director") TraverseCallback;
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}
%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }
}

%include <stl.i>
%include <std_set.i>
%include <std_map.i>
%include <std_shared_ptr.i>;
//%include "std_unique_ptr.i"

%{
#define SWIG_FILE_WITH_INIT
#include "SAScore_data_loader.hpp"
#include "random_seed.hpp"
%}

%include "SAScore_data_loader.hpp"
%include "random_seed.hpp"

// the molpher API
%include "API.i";
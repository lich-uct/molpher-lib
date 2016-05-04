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
%}

%include "SAScore_data_loader.hpp"

// the molpher API
%include "API.i";
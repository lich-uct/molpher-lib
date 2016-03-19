%module(directors="1", threads="1") core // FIXME: not too cool to be releasing GIL on every call (http://swig.10945.n7.nabble.com/How-to-release-Python-GIL-td5027.html)
%feature("director") TreeOperation;
%feature("director") TraverseCallback;
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
#include "SAScore_data_loader.hpp"
%}

%include "SAScore_data_loader.hpp"

// the molpher API
%include "API.i";
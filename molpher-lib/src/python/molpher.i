%module(directors="1", threads="1") core // not top cool to be releasing GIL on every call (http://swig.10945.n7.nabble.com/How-to-release-Python-GIL-td5027.html)
%feature("director") TreeOperation;
%feature("director") TraverseCallback;
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
#include "morphing_functions.hpp"
#include "extensions/SAScore.h"
%}

//%init %{
//    SAScore::loadData();
//%}
%include "../../../backend/extensions/SAScore.h"

// complete morphing function
%include "include/morphing_functions.hpp";

// the molpher API
%include "include/molpher_API/molpher_API.i";
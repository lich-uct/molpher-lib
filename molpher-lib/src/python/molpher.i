%module molpher
%include "std_string.i"

%{
#define SWIG_FILE_WITH_INIT
#include "../../include/morphing_functions.hpp"
%}

void run_path_finder(
    const std::string &storagePath
    , const std::string &jobFile
    , int threadCnt
);

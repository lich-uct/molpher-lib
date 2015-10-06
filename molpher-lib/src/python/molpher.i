%module molpher
%include "std_string.i"

%{
#define SWIG_FILE_WITH_INIT
#include "../../include/testing/testing.hpp"
%}

void run_path_finder(
    const std::string &storagePath
    , const std::string &jobFile
    , int threadCnt
);

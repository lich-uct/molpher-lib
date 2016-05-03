#!/bin/bash

cd $BASE_DIR

LIB_DIR="$PREFIX/lib/"
mkdir -p $LIB_DIR
TBB_LIBS_DIR="${BASE_DIR}/deps/tbb/lib/intel64/gcc4.4/"
cp ${TBB_LIBS_DIR}/libtbbmalloc_proxy.* $LIB_DIR
cp ${TBB_LIBS_DIR}/libtbbmalloc.* $LIB_DIR
cp ${TBB_LIBS_DIR}/libtbb.* $LIB_DIR

INCLUDE_DIR="$PREFIX/include/tbb/"
mkdir -p $INCLUDE_DIR
cp -r deps/tbb/include/tbb/* $INCLUDE_DIR

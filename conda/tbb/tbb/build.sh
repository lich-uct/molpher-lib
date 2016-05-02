#!/bin/bash

cd $BASE_DIR

LIB_DIR="$PREFIX/lib/"
mkdir -p $LIB_DIR
cp dist/lib/libtbbmalloc_proxy.* $LIB_DIR
cp dist/lib/libtbbmalloc.* $LIB_DIR
cp dist/lib/libtbb.* $LIB_DIR

INCLUDE_DIR="$PREFIX/include/tbb/"
mkdir -p $INCLUDE_DIR
cp -r deps/tbb/include/tbb/* $INCLUDE_DIR

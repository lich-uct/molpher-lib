#!/bin/bash

# build and install the library
cd $BASE_DIR
pwd
make -f Makefile CONF=Release_Linux_amd64

LIB_DIR="$PREFIX/lib/"
mkdir -p $LIB_DIR
cp dist/lib/* $LIB_DIR

INCLUDE_DIR="$PREFIX/include/molpher/"
mkdir -p $INCLUDE_DIR
cp -r include/* $INCLUDE_DIR

# build the bindings and install the package
$PYTHON setup.py install

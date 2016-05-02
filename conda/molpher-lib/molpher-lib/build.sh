#!/bin/bash

# build the C++ library
cd $BASE_DIR
pwd
make -f Makefile -j `grep -c ^processor /proc/cpuinfo` CONF=Release_Linux_amd64

# install the C++ library
LIB_DIR="$PREFIX/lib/"
mkdir -p $LIB_DIR
cp dist/lib/libmolpher.* $LIB_DIR

INCLUDE_DIR="$PREFIX/include/molpher/"
mkdir -p $INCLUDE_DIR
rsync -avz --progress include/ $INCLUDE_DIR --exclude '*.i'

# build the bindings and install the Python package
$PYTHON setup.py install

#!/bin/bash

# setup the project
BUILD_DIR=$BASE_DIR/build/molpher-lib/
JOBS=`grep -c ^processor /proc/cpuinfo`
mkdir -p $BUILD_DIR
cd $BUILD_DIR
echo "Building molpher-lib binaries in: `pwd`"
cmake $BASE_DIR -DINSTALL_TBB=OFF

# build and install the C++ library
make -j $JOBS molpher_install

# install the C++ library to the build environment
cp -r $BASE_DIR/dist/. $PREFIX

# build the bindings and install the Python package to the build environment
cd $BASE_DIR
cp res/SAScore.dat src/python/molpher/swig_wrappers/
$PYTHON setup.py install

# clean up
rm -rf ${BASE_DIR}/dist
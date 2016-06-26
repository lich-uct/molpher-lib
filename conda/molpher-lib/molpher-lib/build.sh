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
make -j $JOBS molpher_build_SWIG_Python # needs to be done so that SAScore.dat is included with the package

# install the C++ library to the build environment
cp -r $BASE_DIR/dist/. $PREFIX

# build the bindings and install the Python package to the build environment
cd $BASE_DIR
$PYTHON setup.py install

# clean up
rm -rf ${BASE_DIR}/dist
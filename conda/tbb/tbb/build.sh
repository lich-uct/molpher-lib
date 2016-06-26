#!/bin/bash

# setup the project
BUILD_DIR=$BASE_DIR/build/tbb/
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake $BASE_DIR -DINSTALL_TBB=OFF
make tbb_install

# install to the build environment
cp -r $BASE_DIR/dist/. $PREFIX

# clean up
rm -rf ${BASE_DIR}/dist

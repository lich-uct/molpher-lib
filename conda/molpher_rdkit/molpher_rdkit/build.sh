#!/bin/bash

# setup the project
BUILD_DIR=$BASE_DIR/build/rdkit/
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake $BASE_DIR -DCMAKE_BUILD_TYPE=Release -DINSTALL_TBB=OFF -DINSTALL_RDKit=OFF -DINSTALL_Boost=OFF
make rdkit_install

# install to the build environment
cp -r $BASE_DIR/dist/. $PREFIX

# clean up
rm -rf ${BASE_DIR}/dist

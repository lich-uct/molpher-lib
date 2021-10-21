#!/bin/bash

set -e -x

# setup the project
BUILD_DIR=$BUILD_DIR/cmake/
JOBS=${JOBS:-"$(grep -c ^processor /proc/cpuinfo)"}

mkdir -p $BUILD_DIR
cd $BUILD_DIR
echo "Building molpher-lib package in: `pwd`"
echo "Use maximum of ${JOBS} processes."

# configure
echo "Configuring binary build..."
cmake $BASE_DIR -DCMAKE_BUILD_TYPE=Release -DINCLUDE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX -DINSTALL_TBB=OFF -DINSTALL_RDKit=OFF -DINSTALL_Boost=OFF
echo "Done."

# build and install the C++ library
echo "Building the binaries..."
make -j $JOBS molpher_install
echo "Done."

# build the bindings and install the Python package
echo "Compiling Python wrappers..."
cd $BASE_DIR
"${PYTHON}" setup.py install
echo "Done."

echo "Installation finished."

echo "Cleaning up..."
rm -rvf $BUILD_DIR
echo "All done."
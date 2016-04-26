#!/usr/bin/env bash

# Generate the doxygen XML files for the C++ code
cd doxygen/
rm -rf xml
rm -f *.db
doxygen config.cfg
cd ..

# build the Sphinx documentation
MOLPHER_PATH=../src/python/molpher/
#PYTHONPATH="${MOLPHER_PATH}:${PYTHONPATH}"

sphinx-apidoc -o source/documentation/python/ $MOLPHER_PATH
#make clean
make html

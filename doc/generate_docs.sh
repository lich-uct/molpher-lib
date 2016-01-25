#!/usr/bin/env bash

# Generate the doxygen XML files for the C++ code
cd doxygen/
rm -rf xml
rm -f *.db
doxygen config.cfg
cd ..

# build the Sphinx documentation
sphinx-apidoc -o source/documentation/python/ ../python/molpher/
make clean
make html

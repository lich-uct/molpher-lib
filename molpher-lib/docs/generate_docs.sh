#!/usr/bin/env bash

sphinx-apidoc -f -o source/documentation/python/ ../python/molpher/
cd doxygen/
rm -rf xml
rm -f *.db
doxygen config.cfg
cd ..
make clean
make html
touch source/index.rst
make html

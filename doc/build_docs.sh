#!/usr/bin/env bash

# stop if anything goes wrong
set -e

# Generate the doxygen XML files for the C++ code
cd doxygen/
rm -rf xml
rm -f *.db
doxygen config.cfg
cd ..

# build the Sphinx documentation
make clean
MOLPHER_PATH=../src/python/
sphinx-apidoc -o source/documentation/python/ $MOLPHER_PATH
make html

# upload the docs to GitHub if requested

if [[ "$@" == *"--upload"* ]]
then
    ./upload_docs.sh && echo "Docs uploaded successfully." && exit
    echo && echo "There was an error during upload. See the messages above for more information." 1>&2
fi

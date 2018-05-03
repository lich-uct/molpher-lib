#!/usr/bin/env bash

# stop if anything goes wrong
set -e

# Generate the doxygen XML files for the C++ code
cd doxygen/
rm -rf xml
rm -f *.db
doxygen config.cfg
cd ..

# clean previous build
make clean

# render sample notebooks
MOLPHER_PATH=../src/python/
cd notebooks/
jupyter nbconvert --to html *.ipynb
mv *.html ../source/_static/
cd ..

# configure
sphinx-apidoc -o source/documentation/python/ $MOLPHER_PATH

# make
make html

# upload the docs to GitHub if requested

if [[ "$@" == *"--upload"* ]]
then
    ./upload_docs.sh && echo "Docs uploaded successfully." && exit
    echo && echo "There was an error during upload. See the messages above for more information." 1>&2
fi

# serve docs if required
if [[ "$@" == *"--serve"* ]]
then
    python serve.py
fi

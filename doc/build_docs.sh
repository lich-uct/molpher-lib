#!/usr/bin/bash

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
mkdir -p source/_static
cd notebooks/
jupyter nbconvert --to html *.ipynb
mv *.html ../source/_static/
cd ..

# configure
export PYTHONPATH="../src/python/"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/"
sphinx-apidoc -o source/documentation/python/ ../src/python/

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

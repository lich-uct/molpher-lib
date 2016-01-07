#!/usr/bin/env bash

VERSION=`awk '/^VERSION/{print $NF}' molpher-lib/version.py`
VERSION=`echo "$VERSION" | sed -r 's/[\.]+/_/g'`
VERSION=`echo ${VERSION//[\'\"]/}`
ARCHIVE_NAME="molpher_v$VERSION.tar.gz"
BACKEND="backend"
MOLPHER_LIB="molpher-lib"
DEPENDENCIES="dependencies"
COMMONS="common"
DOCUMENTATION="docs"

# create documentation folder
mkdir -p ${DOCUMENTATION}
cp -r molpher-lib/docs/build/html/ ${DOCUMENTATION}

tar -czvf "../$ARCHIVE_NAME" \
    --exclude=${DEPENDENCIES}/*.sh \
    install-python.sh \
    license.txt \
    ${BACKEND}/auxiliary/ \
    ${BACKEND}/chem/ \
    ${BACKEND}/coord/ \
    ${BACKEND}/core/ \
    ${BACKEND}/dist/ \
    ${BACKEND}/extensions/ \
    ${BACKEND}/nbproject/*.mk \
    ${BACKEND}/TestFiles/ \
    ${BACKEND}/tests/ \
    ${BACKEND}/*.cpp \
    ${BACKEND}/*.h \
    ${BACKEND}/*.hpp \
    ${BACKEND}/Makefile \
    ${BACKEND}/*.dat \
    ${COMMONS}/ \
    ${DEPENDENCIES}/ \
    ${DOCUMENTATION}/ \
    ${MOLPHER_LIB}/lib/ \
    ${MOLPHER_LIB}/include/ \
    ${MOLPHER_LIB}/src/ \
    ${MOLPHER_LIB}/nbproject/*.mk \
    ${MOLPHER_LIB}/python/ \
    ${MOLPHER_LIB}/tests/ \
    ${MOLPHER_LIB}/*.sh \
    ${MOLPHER_LIB}/Makefile \
    ${MOLPHER_LIB}/*.dat \
    ${MOLPHER_LIB}/*.py 2> 'create_linux_release_errors.log'

# cleanup
rm -r ${DOCUMENTATION}
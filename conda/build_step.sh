#!/bin/bash

export BASE_DIR=${BASE_DIR:-"$(pwd)/.."}
export BUILD_DIR=${BUILD_DIR:-"${BASE_DIR}/conda-build"}
mkdir -p $BUILD_DIR
conda build molpher-lib -c conda-forge --python "${PYTHON_VERSION:-3}" --croot "${BUILD_DIR}/conda/"
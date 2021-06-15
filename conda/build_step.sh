#!/bin/bash

export BASE_DIR=${BASE_DIR:-"$(pwd)/.."}
conda build molpher-lib -c conda-forge --python "${PYTHON_VERSION:-3}" --croot "${BASE_DIR}/conda-build/conda/"
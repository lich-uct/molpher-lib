#!/bin/bash

export BASE_DIR="$(pwd)/.."
cp $BASE_DIR/LICENSE.md molpher-lib/
conda build molpher-lib -c conda-forge --python "${PYTHON_VERSION:-3}" --croot "${BASE_DIR}/conda-build/conda/"
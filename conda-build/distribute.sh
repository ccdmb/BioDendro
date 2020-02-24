#!/usr/bin/env bash

set -eu

mkdir -p build

PY_VERSIONS="3.6 3.7 3.8"

for PY_VERSION in ${PY_VERSIONS}
do
    CONDA_OUTPUT=$(conda build --output-folder build --variants "{'python': ['${PY_VERSION}']}" .)
    CREATED_FILE=$(echo "${CONDA_OUTPUT}" | sed -n '/^anaconda upload/s/anaconda upload //p')
    echo "${CREATED_FILE}"
    anaconda upload --all -d "test" --skip-existing "${CREATED_FILE}"
    #conda convert --output-dir "build" --platform all "${CREATED_FILE}"
done

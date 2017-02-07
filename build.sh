#!/usr/bin/env bash

pip install intervaltree sqlalchemy_utils pyfaidx python-magic drmaa snakemake;

minor=$(python -c "import sys; print(sys.version_info.minor)")

${PYTHON} setup.py bdist_wheel || exit 1;

wheel=$(ls dist/*whl);
pip install --no-deps ${wheel}
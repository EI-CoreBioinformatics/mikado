#!/usr/bin/env bash

pip install intervaltree sqlalchemy_utils pyfaidx python-magic drmaa snakemake

minor=$(python -c "import sys; print(sys.version_info.minor)")

python setup.py bdist_wheel
pip install dist/Mikado-1.0.0b10-cp3${minor}-cp3${minor}-*whl;
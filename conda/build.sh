#!/usr/bin/env bash

minor=$(python -c "import sys; print(sys.version_info.minor)")

${PYTHON} setup.py bdist_wheel || exit 1;

wheel=$(ls dist/*3${minor}*whl);
pip install --no-deps ${wheel};

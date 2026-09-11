#! /bin/bash

# propigate linting failure
set -e

# lint for Python syntax errors or undefined names\
ruff check --format=github --select=E9,F63,F7,F82,W605 --extend-exclude=docs/source/schemas/interfaces,src/dendropy/dataio/_beast2_v2_7_8_lark_standalone.py --target-version=py37 .

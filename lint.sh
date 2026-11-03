#! /bin/bash

# propigate linting failure
set -e

# lint for Python syntax errors or undefined names\
ruff check --format=github --select=E9,F63,F7,F82,W605 --extend-exclude=docs/source/schemas/interfaces,src/dendropy/dataio/_beast2_v2_7_8_lark_standalone.py --target-version=py37 .

# F811 reports a name redefined before use, which is what a duplicated test
# method or class name is: Python keeps the last definition and the earlier
# test silently stops running. Scoped to tests/unittests because that tree is
# clean; `ruff check --select=F811 .` currently reports 2 more under
# tests/smoketests and 7 under src/ that need their own decisions.
ruff check --format=github --select=F811 --target-version=py37 tests/unittests

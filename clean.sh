#! /bin/sh

rm -rfv $(find . -name "*.pyc")
rm -rfv $(find . -name "__pycache__")
rm -rfv 'tests/output/'*
rm -rfv 'tests/coverage/'*
rm -rfv src/DendroPy.egg-info
rm -rfv build
rm -rfv dist

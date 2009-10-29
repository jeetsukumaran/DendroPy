#! /bin/sh

rm -rfv $(find . -name "*.pyc")
rm -rfv dendropy/tests/output/*
rm -rfv DendroPy.egg-info
rm -rfv build
rm -rfv dist

#! /bin/sh

python setup.py sdist upload \
    && python setup.py build_sphinx \
    && python setup.py upload_sphinx

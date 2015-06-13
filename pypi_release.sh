#! /bin/sh

python setup.py register sdist upload \
    && python setup.py build_sphinx \
    && python setup.py upload_docs

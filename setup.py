#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Package setup and installation.
"""

import os
import sys
import re
from setuptools import setup, find_packages

def _read(path_components, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *path_components)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get("encoding", "utf8")) as src:
            s = src.read()
        return s

def _read_requirements(path):
    return [
        line.strip()
        for line in _read([path]).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]

project_init = _read(["src", "dendropy", "__init__.py"])
__version__ = re.match(r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)
__project__ = re.match(r".*^__project__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)

setup(
    name=__project__,
    version=__version__,
    author='Jeet Sukumaran, Mark T. Holder, and Matt Moreno',
    author_email='jeetsukumaran@gmail.com, mtholder@ku.edu',
    packages=find_packages("src"),
    package_dir={"": "src"},
    entry_points={
        'console_scripts': [
            # Going forward ...
            'sumtrees=dendropy.application.sumtrees:main',
            'sumlabels=dendropy.application.sumlabels:main',
            'dendropy-format=dendropy.application.dendropy_format:main',
            # Legacy: to be deprecated
            'sumtrees.py=dendropy.application.sumtrees:main',
            'sumlabels.py=dendropy.application.sumlabels:main',
        ],
    },
    include_package_data=True,
    # MANIFEST.in: only used in source distribution packaging.
    # ``package_data``: only used in binary distribution packaging.
    package_data={
        "": [
            "*.txt",
            "*.md",
            "*.rst",
        ],
        "dendropy": [
            # For files in this package's direct namespace
            # (e.g., "src/{normalized_project_name}/*.json")
            # "*.json",
            # For files in a (non-subpackage) subdirectory direct namespace
            # (e.g., "src/{normalized_project_name}/resources/config/*.json")
            # "resources/config/*.json",
            # For files located in 'src/piikun-data/'
            # "../piikun-data/*.json",
            # For files located in 'resources'/'
            # "../../resources/*.json",
        ],
    },
    test_suite = "tests",
    url='http://github.com/jeetsukumaran/DendroPy',
    license='BSD',
    description="A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.",
    long_description=_read(["README.rst"]),
    # long_description_content_type="text/markdown",
    long_description_content_type="text/x-rst",
    # install_requires=_read_requirements("requirements.txt"),
    # extras_require={"test": _read_requirements("requirements-test.txt")},
)


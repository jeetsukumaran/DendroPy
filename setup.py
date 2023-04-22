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

import sys
import re
import os
import io

###############################################################################
### "highly-discouraged" workaround to deal with Python 2.8+ and 3.4+ verifying
### certificates by default; this breaks PyPI uploads ...
import ssl
if hasattr(ssl, '_create_unverified_context'):
    ssl._create_default_https_context = ssl._create_unverified_context

###############################################################################
## Utility

def _read(names, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *names)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get('encoding', 'utf8')) as src:
            s = src.read()
        return s

def _compose_list(values, prefix=None):
    if prefix is None:
        prefix = ""
    s = []
    max_len = max(len(i) for i in values)
    template = "{}{{value:<{}}}".format(prefix,max_len)
    for v in values:
        s.append(template.format(value=v))
    return "\n".join(s)

###############################################################################
# Identification

__version__ = re.match(r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", _read(["src", "dendropy", "__init__.py"]), re.S | re.M).group(1)
sys.stderr.write("-setup.py: DendroPy version {}\n".format(__version__))

###############################################################################
# setuptools/distutils/etc. import and configuration

try:
    import setuptools
    try:
        setuptools_path = " ('" +  os.path.abspath(setuptools.__file__) + "')"
    except OSError:
        setuptools_path = ""
    sys.stderr.write("-setup.py: using setuptools{}\n".format(setuptools_path))
    from setuptools import setup, find_packages
except ImportError as e:
    sys.stderr.write("-setup.py: using distutils\n")
    from distutils.core import setup
    sys.stderr.write("-setup.py: using canned package list\n")
    PACKAGES = [
            "dendropy",
            "dendropy.calculate",
            "dendropy.dataio",
            "dendropy.datamodel",
            "dendropy.interop",
            "dendropy.legacy",
            "dendropy.mathlib",
            "dendropy.model",
            "dendropy.simulate",
            "dendropy.utility",
            "dendropy.utility.libexec",
            ]
else:
    sys.stderr.write("-setup.py: searching for packages\n")
    PACKAGES = find_packages("src")
EXTRA_KWARGS = dict(
    install_requires = ['setuptools'],
    include_package_data = True,
    test_suite = "tests",
    zip_safe = True,
    )

sys.stderr.write("-setup.py: packages identified:\n{}\n".format(_compose_list(PACKAGES, prefix="           - ")))
ENTRY_POINTS = {}

###############################################################################
# Script paths

SCRIPT_SUBPATHS = [
    ['applications', 'sumtrees', 'sumtrees.py'],
    ['applications', 'sumlabels', 'sumlabels.py'],
    ['applications', 'dendropy-format', 'dendropy-format'],
    # ['scripts', 'sumtrees', 'cattrees.py'],
    # ['scripts', 'calculators', 'strict_consensus_merge.py'],
    # ['scripts', 'calculators', 'long_branch_symmdiff.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("-setup.py: scripts identified:\n{}\n".format(_compose_list(SCRIPTS, prefix="           - ")))

###############################################################################
# setuptools/distuils command extensions

try:
    from setuptools import Command
except ImportError:
    sys.stderr.write("-setup.py: setuptools.Command could not be imported: setuptools extensions not available\n")
else:
    sys.stderr.write("-setup.py: setuptools command extensions are available\n")
    command_hook = "distutils.commands"
    ENTRY_POINTS[command_hook] = []

    ###########################################################################
    # coverage
    # from dendropy.test.support import coverage_analysis
    # if coverage_analysis.DENDROPY_COVERAGE_ANALYSIS_AVAILABLE:
    #     sys.stderr.write("-setup.py: coverage analysis available ('python setup.py coverage')\n")
    #     ENTRY_POINTS[command_hook].append("coverage = dendropy.test.support.coverage_analysis:CoverageAnalysis")
    # else:
    #     sys.stderr.write("-setup.py: coverage analysis not available\n")


###############################################################################
# Main setup

### compose long description ###
long_description = _read(["README.rst"])
long_description = long_description.replace("DendroPy-4.x.x", "DendroPy-{}".format(__version__))
long_description = long_description.replace("""download the source code archive""",
    """`download the source code archive <http://pypi.python.org/packages/source/D/DendroPy/DendroPy-{}.tar.gz>`_""".format(__version__))

long_description = long_description + ("""\

Current Release
===============

The current release of DendroPy is version {}.

""".format(__version__))

setup(name='DendroPy',
      version=__version__,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeetsukumaran@gmail.com, mtholder@ku.edu',
      url='http://pypi.org/project/DendroPy//',
      description="A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.",
      license='BSD',
      packages=PACKAGES,
      package_dir={"": "src"},
      # not needed?
      # package_data={
      #     # "dendropy.utility" : ["libexec/*"],
      #     },
      scripts = SCRIPTS,
      long_description=long_description,
      entry_points = ENTRY_POINTS,
      classifiers = [
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics phylogeny phylogenies phylogeography evolution evolutionary biology systematics coalescent population genetics phyloinformatics bioinformatics',
      **EXTRA_KWARGS
      )

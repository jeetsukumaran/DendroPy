#! /usr/bin/env python

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
import os

### "highly-discouraged" workaround to deal with Python 2.8+ and 3.4+ verifying
### certificates by default; this breaks PyPI uploads ...
import ssl
if hasattr(ssl, '_create_unverified_context'):
    ssl._create_default_https_context = ssl._create_unverified_context

###############################################################################
# Identification

from dendropy import __version__, revision_description, description
sys.stderr.write("-setup.py: {}\n".format(description()))

###############################################################################
# setuptools/distutils/etc. import and configuration

try:
    import ez_setup
    try:
        ez_setup_path = " ('" + os.path.abspath(ez_setup.__file__) + "')"
    except OSError:
        ez_setup_path = ""
    sys.stderr.write("-setup.py: using ez_setup{}\n".format(ez_setup_path))
    ez_setup.use_setuptools()
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
            "dendropy.test",
            "dendropy.utility",
            "dendropy.utility.libexec",
            ]
else:
    sys.stderr.write("-setup.py: searching for packages\n")
    PACKAGES = find_packages()
EXTRA_KWARGS = dict(
    install_requires = ['setuptools'],
    include_package_data = True,
    test_suite = "dendropy.test",
    zip_safe = True,
    )

PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
PACKAGE_INFO = [("{p[0]:>40} : {p[1]}".format(p=p)) for p in zip(PACKAGES, PACKAGE_DIRS)]
sys.stderr.write("-setup.py: packages identified:\n{}\n".format("\n".join(PACKAGE_INFO)))
ENTRY_POINTS = {}

###############################################################################
# Script paths

SCRIPT_SUBPATHS = [
    ['applications', 'sumtrees', 'sumtrees.py'],
    # ['scripts', 'sumtrees', 'cattrees.py'],
    # ['scripts', 'sumtrees', 'sumlabels.py'],
    # ['scripts', 'calculators', 'strict_consensus_merge.py'],
    # ['scripts', 'calculators', 'long_branch_symmdiff.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("\n-setup.py: scripts identified: {}\n".format(", ".join(SCRIPTS)))

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
    from dendropy.test.support import coverage_analysis
    if coverage_analysis.DENDROPY_COVERAGE_ANALYSIS_AVAILABLE:
        sys.stderr.write("-setup.py: coverage analysis available ('python setup.py coverage')\n")
        ENTRY_POINTS[command_hook].append("coverage = dendropy.test.support.coverage_analysis:CoverageAnalysis")
    else:
        sys.stderr.write("-setup.py: coverage analysis not available\n")


###############################################################################
# Main setup

### compose long description ###
long_description = open('README.rst').read()
long_description = long_description.replace("DendroPy-4.x.x", "DendroPy-{}".format(__version__))
long_description = long_description.replace("""download the source code archive""",
    """`download the source code archive <http://pypi.python.org/packages/source/D/DendroPy/DendroPy-{}.tar.gz>`_""".format(__version__))

revision_text = revision_description()
long_description = long_description + ("""\

Current Release
===============

The current release of DendroPy is version {}{}.

""".format(__version__, revision_text))

setup(name='DendroPy',
      version=__version__,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeetsukumaran@gmail.com and mtholder@ku.edu',
      url='http://packages.python.org/DendroPy/',
      description="A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.",
      license='BSD',
      packages=PACKAGES,
      package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
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
            "Programming Language :: Python :: 3.1",
            "Programming Language :: Python :: 3.2",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics phylogeny phylogenies phylogeography evolution evolutionary biology systematics coalescent population genetics phyloinformatics bioinformatics',
      **EXTRA_KWARGS
      )

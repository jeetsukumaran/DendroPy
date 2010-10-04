#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
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

###############################################################################
# setuptools/distutils/etc. import and configuration

try:
    import ez_setup
    try:
        ez_setup_path = " ('" + os.path.abspath(ez_setup.__file__) + "')"
    except OSError:
        ez_setup_path = ""
    sys.stderr.write("using ez_setup%s\n" %  ez_setup_path)
    ez_setup.use_setuptools()
    import setuptools
    try:
        setuptools_path = " ('" +  os.path.abspath(setuptools.__file__) + "')"
    except OSError:
        setuptools_path = ""
    sys.stderr.write("using setuptools%s\n" % setuptools_path)
    from setuptools import setup, find_packages
except ImportError, e:
    sys.stderr.write("using distutils\n")
    from distutils.core import setup
    sys.stderr.write("using canned package list\n")
    PACKAGES = ['dendropy',
                'dendropy.dataio',
                'dendropy.dataobject',
                'dendropy.test',
                'dendropy.test.support',
                'dendropy.utility']
    EXTRA_KWARGS = {}
else:
    sys.stderr.write("searching for packages\n")
    PACKAGES = find_packages()
    EXTRA_KWARGS = dict(
        install_requires = ['setuptools'],
        include_package_data=True,
        test_suite = "dendropy.test"
    )

PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
PACKAGE_INFO = [("% 40s : %s" % p) for p in zip(PACKAGES, PACKAGE_DIRS)]
sys.stderr.write("packages identified:\n%s\n" % ("\n".join(PACKAGE_INFO)))
ENTRY_POINTS = {}

###############################################################################
# Script paths

SCRIPT_SUBPATHS = [
    ['scripts', 'sumtrees', 'sumtrees.py'],
    ['scripts', 'sumtrees', 'cattrees.py'],
    ['scripts', 'calculators', 'strict_consensus_merge.py'],
    ['scripts', 'calculators', 'long_branch_symmdiff.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("\nscripts identified: %s\n" % ", ".join(SCRIPTS))

###############################################################################
# setuptools/distuils command extensions

try:
    from setuptools import Command
except ImportError:
    sys.stderr.write("setuptools.Command could not be imported: setuptools extensions not available\n")
else:
    sys.stderr.write("setuptools command extensions are available\n")
    command_hook = "distutils.commands"
    ENTRY_POINTS[command_hook] = []

    ###########################################################################
    # coverage
    from dendropy.test.support import coverage_analysis
    if coverage_analysis.DENDROPY_COVERAGE_ANALYSIS_AVAILABLE:
        sys.stderr.write("coverage analysis available ('python setup.py coverage')\n")
        ENTRY_POINTS[command_hook].append("coverage = dendropy.test.support.coverage_analysis:CoverageAnalysis")
    else:
        sys.stderr.write("coverage analysis not available\n")


###############################################################################
# Main setup

from dendropy import __version__, __revision__
EXTRA_KWARGS["zip_safe"] = True

### compose long description ###
long_description = open('README.txt').read()
long_description = long_description.replace("DendroPy-3.x.x", "DendroPy-%s" % __version__)
long_description = long_description.replace("""download the source code archive""",
    """`download the source code archive <http://pypi.python.org/packages/source/D/DendroPy/DendroPy-%s.tar.gz>`_""" % __version__)

if __revision__.is_available:
    revision_text = " (revision: %s, %s)" % (__revision__.commit_id, str(__revision__.commit_date))
else:
    revision_text = ""
long_description = long_description + ("""\

Current Release
===============

The current release of DendroPy is version %s%s.

""" % (__version__, revision_text))

setup(name='DendroPy',
      version=__version__,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeet@ku.edu and mtholder@ku.edu',
      url='http://packages.python.org/DendroPy/',
      description="A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.",
      license='BSD',
      packages=PACKAGES,
      package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
      # For some reason, following does not work in 2.5 (not tried in 2.6),
      # so this packaging is now implemented through processing of MANIFEST.in
#      package_data={
#        "" : ['doc/Makefile',
#              '/doc/source',
#              'extras'
#             ],
#        "dendropy.test" : ["data/trees"],
#      },
      scripts = SCRIPTS,
      long_description=long_description,
      entry_points = ENTRY_POINTS,
      classifiers = [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics phylogeny phylogenies phylogeography evolution evolutionary biology systematics coalescent population genetics phyloinformatics bioinformatics',
      **EXTRA_KWARGS
      )

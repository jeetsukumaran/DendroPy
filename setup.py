#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Package setup and installation.
"""

import sys
import os

try:
    import ez_setup
    sys.stderr.write("using ez_setup ('%s')\n" % os.path.abspath(ez_setup.__file__))
    ez_setup.use_setuptools()
    import setuptools
    sys.stderr.write("using setuptools ('%s')\n" % os.path.abspath(setuptools.__file__))
    from setuptools import setup, find_packages
except ImportError, e:
    sys.stderr.write("using distutils\n")
    from distutils.core import setup
    PACKAGES = ['dendropy',
                'dendropy.dataio',
                'dendropy.dataobject',
                'dendropy.test',
                'dendropy.test.support',
                'dendropy.utility']
    EXTRA_KWARGS = {}
else:
    PACKAGES = find_packages()
    EXTRA_KWARGS = dict(
        install_requires = ['setuptools'],
        include_package_data=True,
        test_suite = "dendropy.test",
        zip_safe=False,
    )
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
PACKAGE_INFO = [("% 40s : %s" % p) for p in zip(PACKAGES, PACKAGE_DIRS)]
sys.stderr.write("packages identified:\n%s\n" % ("\n".join(PACKAGE_INFO)))

from dendropy import PACKAGE_VERSION

SCRIPT_SUBPATHS = [
    ['scripts', 'sumtrees', 'sumtrees.py'],
    ['scripts', 'sumtrees', 'cattrees.py'],
    ['scripts', 'calculators', 'strict_consensus_merge.py'],
    ['scripts', 'calculators', 'long_branch_symmdiff.py'],
]
SCRIPTS = [os.path.join(*i) for i in SCRIPT_SUBPATHS]
sys.stderr.write("\nscripts: %s\n" % ", ".join(SCRIPTS))

setup(name='DendroPy',
      version=PACKAGE_VERSION,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeet@ku.edu and mtholder@ku.edu',
      url='http://packages.python.org/DendroPy/',
      description="Phylogenetic computing library.",
      license='GPL 3+',
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
      long_description="""\
A Python library for phylogenetic scripting, simulation,
data processing and manipulation.""",
      classifiers = [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics evolution biology',
      **EXTRA_KWARGS
      )

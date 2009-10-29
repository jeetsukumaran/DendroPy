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
                'dendropy.dataobjects',
                'dendropy.tests',
                'dendropy.utility']
    EXTRA_KWARGS = {}
else:
    PACKAGES = find_packages()
    EXTRA_KWARGS = dict(
        install_requires = ['setuptools'],
        include_package_data=True,
        test_suite = "dendropy.tests",
        zip_safe=True,
    )
sys.stderr.write("packages identified:\n    %s\n" % ("\n    ".join(PACKAGES)))
from dendropy import PACKAGE_VERSION
SCRIPT_NAMES = []
setup(name='DendroPy',
      version=PACKAGE_VERSION,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeet@ku.edu and mtholder@ku.edu',
      url='http://packages.python.org/DendroPy/',
      description="Phylogenetic computing library.",
      license='GPL 3+',
      packages=PACKAGES,
      package_data={
        "" : ['doc/*'],
        "dendropy" : ["tests/data/*"]
      },
      scripts = [('scripts/%s' % i) for i in SCRIPT_NAMES],
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

#! /usr/bin/env python

############################################################################
##  setup.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
Package setup and installation.
"""
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup
from setuptools import find_packages
from dendropy import PACKAGE_VERSION

import sys
import os
import subprocess

setup(name='DendroPy',
      version=PACKAGE_VERSION,     
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeet@ku.edu and mtholder@ku.edu',
      url='http://sourceforge.net/projects/dendropy/',
      description="""\
Phylogenetic computing library""",
      license='GPL 3+',
      packages=['dendropy'],
      package_dir={'dendropy': 'dendropy'},
      package_data={
        "" : ['doc/*'],
        "dendropy" : ["tests/data/*"]
      },
      scripts=['scripts/sumtrees.py'],   
      test_suite = "dendropy.tests",
      include_package_data=True,         
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,      
      long_description="""\
A Python library for phylogenetic computation, scripting, simulation,
data processing and manipulation.""",
      classifiers = [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU Library or  General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics evolution biology',      
      )

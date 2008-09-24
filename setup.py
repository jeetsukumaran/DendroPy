#! /usr/bin/env python

############################################################################
##  setup.py
##
##  Part of the DendroPy phylogenetic library.
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Package setup and installation.
"""

from setuptools import setup, find_packages
import sys, os

version = '2.0.0'

setup(name='DendroPy',
      version=version,
      author='Jeet Sukumaran and Mark T. Holder',
      author_email='jeetsukumaran@gmail.com and mtholder@ku.edu',
      maintainer = "Jeet Sukumaran and Mark Holder", 
      maintainer_email = "jeetsukumaran@frogweb.org mtholder@ku.edu",      
      description="Phylogenetic library",
      url='',
      license='GPL 3+',
      packages=find_packages(), #find_packages(exclude=['ez_setup', 'examples', 'tests']),
      package_data={
        "dendropy" : ["tests/sources/*"]
      },
      include_package_data=True,      
      test_suite = "dendropy.tests",
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,      
      long_description="""\
A pure-Python library for the serialization/deserialization, parsing,
analysis, manipulations and other operations on phylogenetic trees and
data.""",
      classifiers = [
            "Development Status :: 2 - Pre-Alpha",
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU Library or  General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
      keywords='phylogenetics bioinformatics tree',      
      )

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
DendroPy testing suite.
"""

import unittest
import re
import os

def get_test_file_names():
    """Get list of test file names."""
    path = os.path.dirname(__file__)
    files = os.listdir(path)
    test_file_pattern = re.compile("test.*\.py$", re.IGNORECASE)
    test_files = []
    for f in files:
        if test_file_pattern.search(f):
            test_files.append("dendropy.test." + os.path.splitext(f)[0])
    return test_files

def get_test_suite(test_file_names=None):
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.test`. Right now, assumes (a) no subdirectories (though
    this can easily be accommodated) and (b) every test to be run is
    sitting in a module with a file name of 'test*.py', and, conversely,
    every file with a name of 'test*.py' has test(s) to be run.
    """
    if test_file_names is None:
        test_file_names = get_test_file_names()
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_file_names)
    return unittest.TestSuite(tests)

def run():
    "Runs all of the unittests"
    runner = unittest.TextTestRunner()
    runner.run(get_test_suite())

if __name__ == "__main__":
    run()


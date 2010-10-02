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


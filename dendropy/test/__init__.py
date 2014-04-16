#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
from distutils.util import strtobool
from dendropy.utility import metavar

def discover_test_module_paths(filter_patterns=None):
    """
    Discovers test modules. If `filter_patterns` is `None`, then
    all files in *immediate* directory that begin with 'test' will
    be added to the set returned. If `filter_patterns` is not `None`, then it
    should be a list of regular expression patterns, and only files that match
    at least one of the patterns will be returned.
    """
    test_module_pattern = re.compile("^test.*\.py$", re.IGNORECASE)
    if filter_patterns:
        filter_pattern = re.compile("(" + r"\|".join(filter_patterns) + ")")
    else:
        filter_pattern = None
    path = os.path.dirname(__file__)
    filenames = os.listdir(path)
    test_modules = []
    for filename in filenames:
        if test_module_pattern.match(filename):
            if filter_pattern is None or filter_pattern.match(filename):
                test_modules.append("dendropy.test." + os.path.splitext(filename)[0])
    return test_modules

def get_test_suite(test_names=None):
    """
    If `test_names` is not `None`, creates a test suite out of those
    modules. Otherwise, creates a test suite from all of the modules in
    `dendropy.test` using the discovery.
    """
    if test_names is None:
        test_names = discover_test_module_paths()
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_names)
    return unittest.TestSuite(tests)

def run(test_suite=None, verbosity=1, failfast=False):
    "Runs all of the unittests"
    runner = unittest.TextTestRunner(verbosity=verbosity, failfast=failfast)
    if test_suite is None:
        test_suite = get_test_suite()
    runner.run(test_suite)

def fail_incomplete_tests():
    return bool(strtobool(os.environ.get(metavar.FAIL_INCOMPLETE_TESTS_ENVAR, "0")))

if __name__ == "__main__":
    run()


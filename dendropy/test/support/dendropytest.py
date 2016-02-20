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
Extension to the basic unittest TestCase.
"""

import collections
import sys
import re
import os
import unittest
from distutils.util import strtobool
from dendropy.utility import metavar
from dendropy.utility import messaging

# Defining this here means that unittest will exclude all lines from this
# module in the traceback report when an assertion fails, allowing
# for the starting point of the traceback to be the point where the assertion
# was made, rather than the point where an exception was raised because
# the assertion was false.
__unittest = True

def discover_test_module_paths(filter_patterns=None):
    """
    Discovers test modules. If ``filter_patterns`` is |None|, then
    all files in *immediate* directory that begin with 'test' will
    be added to the set returned. If ``filter_patterns`` is not |None|, then it
    should be a list of regular expression patterns, and only files that match
    at least one of the patterns will be returned.
    """
    test_module_pattern = re.compile("^test.*\.py$", re.IGNORECASE)
    if filter_patterns:
        filter_pattern = re.compile("(" + r"\|".join(filter_patterns) + ")")
    else:
        filter_pattern = None
    path = os.path.dirname(os.path.dirname(__file__))
    filenames = os.listdir(path)
    test_modules = []
    for filename in filenames:
        if test_module_pattern.match(filename):
            if filter_pattern is None or filter_pattern.match(filename):
                # test_modules.append("" + os.path.splitext(filename)[0])
                test_modules.append("dendropy.test." + os.path.splitext(filename)[0])
    return test_modules

def get_test_suite(test_names=None):
    """
    If ``test_names`` is not |None|, creates a test suite out of those
    modules. Otherwise, creates a test suite from all of the modules in
    ``dendropy.test`` using the discovery.
    """
    if test_names is None:
        test_names = discover_test_module_paths()
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_names)
    return unittest.TestSuite(tests)

class ExtendedTestCase(unittest.TestCase):
    """
    Extends unittest.TestCase with various new assertion tests.
    """

    def _get_logger(self):
        if not hasattr(self, "_logger") or self._logger is None:
            self._logger = messaging.get_logger(self.__class__.__name__)
        return self._logger
    def _set_logger(self, logger):
        self._logger = logger
    logger = property(_get_logger, _set_logger)

    def assertCountEqual(self, *args, **kwargs):
        if sys.hexversion >= 0x03020000:
            super(ExtendedTestCase, self).assertCountEqual(*args, **kwargs)
        else:
            self.assertEqual(collections.Counter(args[0]), collections.Counter(args[1]))

    def fail_incomplete_tests(self):
        return bool(strtobool(os.environ.get(metavar.FAIL_INCOMPLETE_TESTS_ENVAR, "0")))

    def assertEqualUnorderedSequences(self, x1, x2):
        c1 = collections.Counter(x1)
        c2 = collections.Counter(x2)
        return self.assertEqual(c1, c2)

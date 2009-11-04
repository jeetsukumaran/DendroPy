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
import os
import re
import sys

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

from dendropy.utility import texttools

###############################################################################
## FILE PATH MAPPING

try:
    import pkg_resources
    _LOG.info("Using local pkg_resources path mapping")
    TESTS_DIR = pkg_resources.resource_filename("dendropy", "tests")
    SCRIPTS_DIR = pkg_resources.resource_filename("dendropy", "scripts")
except:
    _LOG.info("Using local filesystem path mapping")
    TESTS_DIR = os.path.dirname(__file__)
    SCRIPTS_DIR = os.path.join(os.path.pardir, "scripts")

TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")
TESTS_OUTPUT_DIR = os.path.join(TESTS_DIR, "output")

###############################################################################
## TESTING LEVELS

class TestLevel:
    FAST, NORMAL, SLOW, EXHAUSTIVE = 0, 10, 20, 30
    def name(i):
        if i <= TestLevel.FAST:
            return "FAST"
        if i <= TestLevel.NORMAL:
            return "NORMAL"
        if i <= TestLevel.SLOW:
            return "SLOW"
        return "EXHAUSTIVE"
    name = staticmethod(name)
    def name_to_int(l):
        try:
            return int(l)
        except:
            pass
        l = l.upper()
        if l == "FAST":
            return TestLevel.FAST
        if l == "NORMAL":
            return TestLevel.NORMAL
        if l == "SLOW":
            return TestLevel.SLOW
        if l == "EXHAUSTIVE":
            return TestLevel.EXHAUSTIVE
        raise ValueError("TestLevel %s unrecognized" % l)
    name_to_int = staticmethod(name_to_int)

def fast_testing_notification(logger, module_name, message=None, level=TestLevel.FAST):
    if message is None:
        message = "tests skipped"
    logger.warning('\nRunning in %s Testing Level. Skipping %s tests in %s: %s' \
        % (TestLevel.name(get_current_testing_level()),
           TestLevel.name(level),
           module_name,
           message))

def get_current_testing_level():
    l = os.environ.get("DENDROPY_TESTING_LEVEL")
    if l is None:
        if "DENDROPY_FAST_TESTS" in os.environ:
            return TestLevel.FAST
        return TestLevel.NORMAL
    try:
        return TestLevel.name_to_int(l)
    except:
        _LOG.warn("the value %s for DENDROPY_TESTING_LEVEL is not recognized.  Using NORMAL level" % l)
    return TestLevel.NORMAL

def is_test_enabled(level, logger=None, module_name="", message=None):
    tl = get_current_testing_level()
    if level > tl:
        if logger:
            fast_testing_notification(logger, module_name, message, level)
        return False
    return True

###############################################################################
## FILE PATHS

def data_source_stream(filename):
    return open(data_source_path(filename), "r")

def tree_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(TESTS_DATA_DIR, "trees", filename)

def data_source_path(filename=None):
    if filename is None:
        filename = ""
    elif isinstance(filename, list):
        filename = os.path.sep.join(filename)
    return os.path.join(TESTS_DATA_DIR, filename)

def named_output_stream(filename=None, suffix_timestamp=True):
    return open(named_output_path(filename=filename, suffix_timestamp=suffix_timestamp), "w")

def named_output_path(filename=None, suffix_timestamp=True):
    if filename is None:
        filename = ""
    else:
        if isinstance(filename, list):
            filename = os.path.sep.join(filename)
        if suffix_timestamp:
            filename = "%s.%s" % (filename, texttools.pretty_timestamp())
    if not os.path.exists(TESTS_OUTPUT_DIR):
        os.makedirs(TESTS_OUTPUT_DIR)
    return os.path.join(TESTS_OUTPUT_DIR, filename)

def scripts_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(SCRIPTS_DIR, filename)

###############################################################################
## TEST SUITES

def get_test_suite():
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.tests`. Right now, assumes (a) no subdirectories (though
    this can easily be accommodated) and (b) every test to be run is
    sitting in a module with a file name of 'test*.py', and, conversely,
    every file with a name of 'test*.py' has test(s) to be run.
    """
    # get list of test file names'
    path = os.path.dirname(__file__)
    files = os.listdir(path)
    test_file_pattern = re.compile("test.*\.py$", re.IGNORECASE)
    test_files = []
    for f in files:
        if test_file_pattern.search(f):
            test_files.append("dendropy.tests." + os.path.splitext(f)[0])

    # extract the tests
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_files)

    # return the suite
    return unittest.TestSuite(tests)

def run():
    "Runs all of the unittests"
    runner = unittest.TextTestRunner()
    runner.run(get_test_suite())

if __name__ == "__main__":
    unittest.main()


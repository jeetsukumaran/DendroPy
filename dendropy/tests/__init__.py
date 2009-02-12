#! /usr/bin/env python

############################################################################
##  __init__.py
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
dendropy testing suite
"""

import unittest
import os
import re
import sys

from dendropy import get_logger
from dendropy.utils import find_files
_LOG = get_logger("tests")
    

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
    logger.warning('\nRunning in %s Testing Level. Skipping %s tests in %s: %s' % (TestLevel.name(get_current_testing_level()), TestLevel.name(level), module_name, message))

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

def data_source_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(os.path.join(os.path.dirname(__file__),'data'), filename)
    
def scripts_source_path(filename=None):
    if filename is None:
        filename = ""
    td = os.path.dirname(__file__)
    dd = os.path.dirname(td)
    return os.path.join(os.path.join(dd,'scripts'), filename)
    
def data_source_trees(format="*", heavy=False):
    if heavy:
        path = data_source_path(heavy)
    else:
        path = data_source_path()
    filename_filter = "*.trees." + format
    files = utils.find_files(top=path,
                                recursive=False,
                                filename_filter=filename_filter,
                                dirname_filter=None,
                                excludes=None,
                                complement=False,
                                respect_case=False,
                                expand_vars=True,
                                include_hidden=False)
    return files                                


def data_target_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(os.path.join(os.path.dirname(__file__),'output'), filename)
    

###############################################################################    
# this is/was a major PIA ...

from decimal import Decimal, ROUND_HALF_UP

def symmetric_arithmetic_round(n, digits=0):
    """
    Symmetric Arithmetic Rounding for decimal numbers
    
    Modified from:
        http://pyxx.org/2007/10/28/how-to-round-decimal-numbers-in-python/
    
    d       - Decimal number to round
    digits  - number of digits after the point to leave
    """
    d = Decimal("%s" % n)
    dval = d.quantize(Decimal("1") / (Decimal('10') ** digits), \
        ROUND_HALF_UP) 
        
def is_almost_equal(n1, n2, precision=2):
    return symmetric_arithmetic_round(n1, precision) \
        == symmetric_arithmetic_round(n2, precision)
        
###############################################################################        
    

def get_test_suite():
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.tests`. Right now, assumes (a) no subdirectories (though
    this can easily be accomodated) and (b) every test to be run is
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


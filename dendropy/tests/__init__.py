#! /usr/bin/env python

############################################################################
##  __init__.py
##
##  Part of the DendroPy phylogenetic computation library.
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
dendropy testing suite
"""

import unittest
import os

from dendropy import get_logger

          
# pylint: disable-msg=C0111,W0401,W0611

def get_test_messenger(name):
    logger = get_logger(name)
    return _LOG.info
    
def test_data_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(os.path.join(os.path.dirname(__file__),'data'), filename)
    
def test_target_path(filename=None):
    if filename is None:
        filename = ""
    return os.path.join(os.path.join(os.path.dirname(__file__),'output'), filename)

def even_more_tests(all_suites):
    "finds test from module introspection"
    if __name__ != "__main__":
         return 
    #commented out
    for i in __all__:
        module = __import__("dendropy.tests.%s" % i)
        _LOG.debug(i)
        tests_mod = getattr(module, "tests")
        sub_test_mod = getattr(tests_mod, i)
        suite = sub_test_mod.additional_tests()
        if suite:
            all_suites.append(suite)

if __name__ == "__main__":
    test_all()


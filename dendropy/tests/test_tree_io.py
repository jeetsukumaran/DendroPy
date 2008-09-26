#! /usr/bin/env python

############################################################################
##  test_tree_io.py
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
Tests input/output of trees from files.
"""

import unittest
import datetime
import logging
import os

from dendropy import get_logger
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("test_tree_io")

from dendropy import utils

### MODULE THAT WE ARE TESTING ###
from dendropy import dataio
### MODULE THAT WE ARE TESTING ###

def iterate_on_trees(tree_files, tf_iterator=dataio.iterate_over_trees):
    logging_level = get_logging_level()
    total_tree_files = len(tree_files)
    total_trees = 0    
    start_time = datetime.datetime.now()
    if logging_level > logging.INFO:
        minimal_logging = True
    else:
        minimal_logging = False
    for tree_file_idx, tree_file in enumerate(tree_files):
        if not minimal_logging:
            _LOG.info("\n  Iterator: %s" % tf_iterator.__name__)
            _LOG.info("      File: %s" % os.path.basename(tree_file))
        for tree_idx, tree in enumerate(tf_iterator(filepath=tree_file)):
            if not minimal_logging:
                _LOG.debug("\n%s" % str(tree))
        total_trees += tree_idx
    if not minimal_logging:        
        _LOG.debug("\n")
    end_time = datetime.datetime.now() 
    _LOG.info("Trees Read: %s" % total_trees)        
    _LOG.info("Start time: %s" % start_time)        
    _LOG.info("  End time: %s" % end_time)
    run_time = end_time-start_time
    _LOG.info("  Run time: %s" % utils.pretty_print_deltatime(run_time))
    return run_time
        
def compare_parse_methods(tree_files, methods):
    results = {}
    for method in methods:
        results[method] = iterate_on_trees(tree_files=tree_files, tf_iterator=method)
    _LOG.info("\n--- TREE ITERATION PERFORMANCE COMPARISON ---")        
    for m1 in methods:
        for m2 in methods[methods.index(m1)+1:]:
            _LOG.info("<%s> vs. <%s> = %s " % (m1.__name__, m2.__name__, results[m1]-results[m2])) 

    
class TreeIOTest(unittest.TestCase):

    def test_newick(self):
        sources = utils.find_files(top=dendropy.tests.test_source_path(),
                                    recursive=False,
                                    filename_filter="*.newick.tre",
                                    dirname_filter=None,
                                    excludes=None,
                                    complement=False,
                                    respect_case=False,
                                    expand_vars=True,
                                    include_hidden=False)
        compare_parse_methods(sources, [dataio.tree_iter, dataio.iterate_over_trees])
        


def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(TreeIOTest)

def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
    """
    return additional_tests()


if __name__ == "__main__":
    unittest.main()
    
    

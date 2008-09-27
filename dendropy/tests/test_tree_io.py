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
from optparse import OptionGroup
from optparse import OptionParser

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
    _LOG.info("\n*** ITERATOR: <%s> ***" % tf_iterator.__name__)
    for tree_file_idx, tree_file in enumerate(tree_files):
        _LOG.info("   - %s" % os.path.basename(tree_file))
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
    _LOG.info("  Run time: %s" % utils.pretty_print_timedelta(run_time))
    return run_time

def compare_parse_performance(tree_files, methods):
    _LOG.info("\nRunning iterators for (speed) performance comparison ...")
    results = {}
    for method in methods:
        results[method] = iterate_on_trees(tree_files=tree_files, tf_iterator=method)
    _LOG.info("\n---")
    for m1 in methods:
        for m2 in methods[methods.index(m1)+1:]:
            t1 = results[m1]
            t2 = results[m2]
            if t1 >= t2:
                diff = t1 - t2
                diff_sign = "+"
            else:
                diff = t2 - t1
                diff_sign = "-"
            diff_seconds = diff.seconds + float(diff.microseconds)/1000000
            _LOG.info("<%s> vs. <%s> = %s%s seconds " % (m1.__name__, m2.__name__, diff_sign, diff_seconds))

def test_tree_iter_performance(filename_filter,
                               heavy=False,
                               wait_to_start=False):
    if heavy:
        top = dendropy.tests.test_data_path('heavy')
    else:
        top = dendropy.tests.test_data_path()
    sources = utils.find_files(top=top,
                                recursive=False,
                                filename_filter=filename_filter,
                                dirname_filter=None,
                                excludes=None,
                                complement=False,
                                respect_case=False,
                                expand_vars=True,
                                include_hidden=False)
    if wait_to_start:
        raw_input("Hit [ENTER] to begin iterating over trees: ")

    iterate_on_trees(sources)

class TreeIOTest(unittest.TestCase):

    def test_dummy(self):
        _LOG.warning("\n\n*** TODO: Check correctness of tree!\n")

def main_local():
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True)

    parser.add_option('-p', '--performance',
                      action="store_const",
                      dest="format_type",
                      const="all",
                      help="evaluate NEXUS and NEWICK parsing performance")

    parser.add_option('--NEXUS',
                      action="store_const",
                      dest="format_type",
                      const="NEXUS",
                      help="evaluate NEXUS format parsing performance")

    parser.add_option('--NEWICK',
                      action="store_const",
                      dest="format_type",
                      const="NEWICK",
                      help="evaluate NEWICK format parsing performance")

    parser.add_option('-H', '--heavy',
                        action='store_true',
                        dest='heavy',
                        default=False,
                        help='run heavy (large file) version of performance tests')

    parser.add_option('-w', '--wait',
                        action='store_true',
                        dest='wait',
                        default=False,
                        help='wait for user confirmation before starting runs')

    (opts, args) = parser.parse_args()

    if opts.format_type == "NEXUS":
        filename_filter = "*.nex.tre"
    elif opts.format_type == "NEWICK":
        filename_filter = "*.newick.tre"
    else:
        filename_filter = "*.tre"
    test_tree_iter_performance(filename_filter=filename_filter,
                               heavy=opts.heavy,
                               wait_to_start=opts.wait)


import sys
if __name__ == "__main__":
    if len(sys.argv) == 1:
        unittest.main()
    else:
        main_local()


    #compare_heavy(dataio.iterate_over_trees, "*.newick.tre")
    #compare_heavy(dataio.tree_iter, "*.newick.tre")

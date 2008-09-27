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
import tempfile
import os
from optparse import OptionGroup
from optparse import OptionParser

from dendropy import get_logger
from dendropy import get_logging_level

import dendropy.tests
_LOG = get_logger("test_tree_io")

from dendropy import taxa
from dendropy import trees
from dendropy import utils

### MODULES THAT WE ARE TESTING ###
from dendropy import dataio
from dendropy import nexus
### MODULES THAT WE ARE TESTING ###

def iterate_on_trees(tree_files, tf_iterator=dataio.iterate_over_trees):
    """
    Test (supposedly) memory-economical iteration on trees.
    """
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
    """
    Test speed of (supposedly) memory-economical iteration on trees.
    """                               
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
    
def get_anolis_consensus_tree():
    leaves = {
            "Anolis_ahli": 0.2642,
            "Anolis_aliniger": 0.16,
            "Anolis_alutaceus": 0.1619,
            "Anolis_angusticeps": 0.0857,
            "Anolis_bahorucoensis": 0.2267,
            "Anolis_barahonae": 0.2115,
            "Anolis_brevirostris": 0.1801,
            "Anolis_coelestinus": 0.1932,
            "Anolis_cristatellus": 0.2144,
            "Anolis_cuvieri": 0.1687,
            "Anolis_distichus": 0.1151,
            "Anolis_equestris": 0.0227,
            "Anolis_garmani": 0.1068,
            "Anolis_grahami": 0.0864,
            "Anolis_insolitus": 0.2439,
            "Anolis_krugi": 0.1573,
            "Anolis_lineatopus": 0.1957,
            "Anolis_loysiana": 0.1836,
            "Anolis_luteogularis": 0.0306,
            "Anolis_marcanoi": 0.2359,
            "Anolis_occultus": 0.4231,
            "Anolis_olssoni": 0.2569,
            "Anolis_ophiolepis": 0.0945,
            "Anolis_paternus": 0.0595,
            "Anolis_sagrei": 0.0968,
            "Anolis_strahmi": 0.1978,
            "Anolis_stratulus": 0.1973,
            "Anolis_valencienni": 0.1643,
            "Anolis_vanidicus": 0.206,
            "Diplolaemus_darwinii": 0.3182,
    }
    
    taxa_block = taxa.TaxaBlock()
    leaf_nodes = []
    for tax_label in leaves:
        taxon = taxa_block.add_taxon(elem_id="TAXON_"+tax_label, label=tax_label)
        node = trees.Node(elem_id="tax_label_"+tax_label, taxon=taxon)
        node.edge.length = leaves[tax_label]

def read_newick_tree(tree_filepath):
    """
    Wrapper to read and return a tree from a single-tree NEWICK file.
    """
    f = open(tree_filepath, 'r')
    tstr = f.read()
    _LOG.info('Reading "%s"' % os.path.basename(tree_filepath))
    _LOG.debug(tstr)
    tree = nexus.parse_newick_string(tstr)
    leaves = tree.leaves()
    _LOG.info("%d leaves on tree: %s" % (len(leaves), (", ".join([str(n.taxon) for n in leaves]))))
    return tree
    
def write_newick_tree(tree, tree_filepath):
    """
    Wrapper to write a single tree to a NEWICK file.
    """
    nw = nexus.NewickTreeWriter()    
    f = open(tree_filepath, 'w')
    tstr = nw.compose_tree(tree)
    _LOG.info('\nWriting "%s"' % os.path.basename(tree_filepath))
    _LOG.debug(tstr)
    f.write(tstr)

class TreeIOTest(unittest.TestCase):

    def parse_tree_file(self,
                        tree_filepath, 
                        reader_method, 
                        writer_method):
        """
        Reads a (single) tree from a (single-)tree file,
        writes it out, and reads it back in again.
        """
        tree1 = reader_method(tree_filepath)
        taxa_block = tree1.infer_taxa_block()
        tree1_fpath = tempfile.NamedTemporaryFile().name
        writer_method(tree1, tree1_fpath)
        tree2 = reader_method(tree1_fpath)
        tree2_fpath = tempfile.NamedTemporaryFile().name
        writer_method(tree2, tree2_fpath)
        tree3 = reader_method(tree2_fpath)
        
        s1 = open(tree1_fpath, "r").read()
        s2 = open(tree1_fpath, "r").read()
        self.assertEqual(s1, s2, "Reparsed tree strings do not match:\n\n%s\n\n%s" % (s1, s2)) 
        _LOG.info("\nReparsed tree string match.")
        
    def test_tree_file_parse(self):    
        self.parse_tree_file(dendropy.tests.test_data_path('anolis.mbcon.newick.tre'),
                             read_newick_tree,
                             write_newick_tree)

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

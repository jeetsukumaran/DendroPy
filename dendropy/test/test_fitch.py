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
Tests of Fitch algorithm calculations.
"""

import random
import unittest
import math
import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

import dendropy
from dendropy.calculate.treescore import fitch_down_pass
from dendropy.test.support import pathmap

class FitchTest(unittest.TestCase):

    def test_pscores_with_gaps_as_new_state(self):
        # #NEXUS
        # begin paup;
        #     set warnroot = no;
        #     exe apternodus.chars.nexus;
        #     gett file = apternodus.tre;
        #     set criterion = parsimony;
        #     pset gap = newstate;
        #     pscore;
        # end;
        expected_scores = [396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 396, 713, 715, 723, 733, 672, 719, 734, 709, 695, 686]
        self.verify_pscores("apternodus.chars.nexus", "apternodus.tre", False, expected_scores)

    def test_pscores_with_gaps_as_missing(self):
        # #NEXUS
        # begin paup;
        #     set warnroot = no;
        #     exe apternodus.chars.nexus;
        #     gett file = apternodus.tre;
        #     set criterion = parsimony;
        #     pset gap = missing;
        #     pscore;
        # end;
        expected_scores = [ 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 671, 670, 678, 687, 633, 675, 689, 668, 652, 644]
        self.verify_pscores("apternodus.chars.nexus", "apternodus.tre", True, expected_scores)


    def verify_pscores(self, char_fname, trees_fname, gaps_as_missing, expected_scores):
        dataset = dendropy.DataSet.get_from_path(
                pathmap.char_source_path(char_fname),
                "nexus")
        dataset.read_from_path(
                pathmap.tree_source_path(trees_fname),
                schema='NEXUS',
                taxon_namespace=dataset.taxon_namespaces[0])
        char_mat = dataset.char_matrices[0]
        # sa = char_mat.default_state_alphabet
        # for x in sa:
        #     print("{}: {}".format(x, x.is_gap_state))
        # for x in sa:
        #     print("{}\t{}\t{}\t\t\t\t{}".format(x, x._index, x.fundamental_indexes, x.fundamental_indexes_with_gaps_as_missing))
        taxon_state_sets_map = char_mat.taxon_state_sets_map(gaps_as_missing=gaps_as_missing)
        tree_list = dataset.tree_lists[0]
        self.assertEqual(len(expected_scores), len(tree_list))
        for n, tree in enumerate(tree_list):
            node_list = tree.postorder_node_iter()
            pscore = fitch_down_pass(node_list, taxon_state_sets_map=taxon_state_sets_map)
            # print("{} vs. {}".format(expected_scores[n], pscore))
            self.assertEqual(expected_scores[n], pscore)

if __name__ == "__main__":
    unittest.main()


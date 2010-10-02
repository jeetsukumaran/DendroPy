
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
Tests of tree metrics.
"""

import random
import unittest
import math
from cStringIO import StringIO

import dendropy
from dendropy.treecalc import fitch_down_pass
from dendropy.test.support import pathmap

from dendropy import treecalc

class FitchTest(unittest.TestCase):

    def testPScore(self):
        expected_scores = [370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 671, 670, 678, 687, 633, 675, 689, 668, 652, 644]
        dataset = dendropy.DataSet(stream=open(pathmap.char_source_path("apternodus.chars.nexus"), "rU"),
                                   schema='NEXUS')
        dataset.read(stream=open(pathmap.tree_source_path("apternodus.tre"), "rU"),
                     schema='NEXUS',
                     taxon_set=dataset.taxon_sets[0])
        char_mat = dataset.char_matrices[0]
        taxa_to_state_set_map = char_mat.create_taxon_to_state_set_map()
        tree_list = dataset.tree_lists[0]
        self.assertEqual(len(expected_scores), len(tree_list))
        for n, tree in enumerate(tree_list):
            node_list = tree.postorder_node_iter()
            pscore = fitch_down_pass(node_list, taxa_to_state_set_map=taxa_to_state_set_map)
            self.assertEqual(expected_scores[n], pscore)

if __name__ == "__main__":
    unittest.main()


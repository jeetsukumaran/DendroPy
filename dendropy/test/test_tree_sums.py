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
Tests of summarization.
"""

import unittest
import dendropy
from dendropy import treecalc
from dendropy import treesum
from dendropy import treesplit
from dendropy.test.support import pathmap

class TestConsensusTree(unittest.TestCase):

    def setUp(self):
        self.tree_list = dendropy.TreeList()
        for t in xrange(1, 5):
            tf = pathmap.tree_source_path('pythonidae.mb.run%d.t' % t)
            self.tree_list.read_from_path(tf, 'nexus', tree_offset=25)
        self.mb_con_tree = dendropy.Tree.get_from_path(
                pathmap.tree_source_path("pythonidae.mb.con"),
                schema="nexus",
                index=0,
                taxon_set=self.tree_list.taxon_set)
        self.mb_con_tree.update_splits()

    def testConsensus(self):
        con_tree = self.tree_list.consensus(min_freq=0.50, trees_splits_encoded=False, support_label_decimals=2)
        con_tree.update_splits()
        self.assertEqual(treecalc.symmetric_difference(self.mb_con_tree, con_tree), 0)
        self.assertEqual(len(con_tree.split_edges), len(self.mb_con_tree.split_edges))
        sd = self.tree_list.split_distribution
        for split in self.mb_con_tree.split_edges:
            edge1 = self.mb_con_tree.split_edges[split]
            edge2 = con_tree.split_edges[split]
            if edge1.head_node.label and edge2.head_node.label:
                s1 = float(edge1.head_node.label)
                s2 = round(float(edge2.head_node.label), 2)
                self.assertAlmostEqual(s1, s2, 2)

if __name__ == "__main__":
    unittest.main()

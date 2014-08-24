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
Tests native tree structuring routines.
"""

import unittest
from dendropy.test.support import pathmap
from dendropy.utility import messaging
from dendropy.test.support.dendropytest import ExtendedTestCase
import dendropy
import re

_LOG = messaging.get_logger(__name__)

class TruncateTree(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path("pythonidae.reference-trees.nexus"), "nexus")

    def check_ultrametric_tree(self, tree, dist):
        self.assertTrue(tree._debug_tree_is_valid())
        tree.calc_node_root_distances()
        for nd in tree.leaf_node_iter():
            self.assertAlmostEqual(nd.root_distance, nd.distance_from_root())
            self.assertAlmostEqual(dist, nd.root_distance)

    def test_truncate_ultrametric(self):
        for tree in self.trees:
            dists = tree.calc_node_root_distances()
            min_dist, max_dist = tree.minmax_leaf_distance_from_root()
            trunc_dists = [(max_dist * f) for f in (0.25, 0.5, 0.75, 0.90)]
            for td in trunc_dists:
                working = dendropy.Tree(tree)
                working.truncate_from_root(td)
                for idx, leaf in enumerate(working.leaf_node_iter()):
                    if leaf.label is None and leaf.taxon is None:
                        leaf.taxon = dendropy.Taxon(label="t%s" % (idx+1))
                self.check_ultrametric_tree(working, td)

class TestTreeLadderization(unittest.TestCase):

    def setUp(self):
        self.tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

    def clean_newick_str(self, s):
        """
        Strips out everything but the core tree statement characters from a
        NEWICK string.
        """
        return re.sub(r'[^()A-H,]', '', s)

    def testLadderizeLeft(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=True)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
                self.clean_newick_str("[&R] ((F,(G,H)),(A,(B,(C,(D,E)))));"))

    def testLadderizeRight(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=False)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
               self.clean_newick_str("[&R] (((((D,E),C),B),A),((G,H),F));"))

class TreeMidpointRootingTest(ExtendedTestCase):

    def testMidpointRooting(self):
        taxa = dendropy.TaxonNamespace()
        test_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.randomly-rooted.tre'),
                "nexus",
                taxon_namespace=taxa,
                rooting="force-rooted")
        expected_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.midpoint-rooted.tre'),
                "nexus",
                taxon_namespace=taxa,
                rooting="force-rooted")
        for idx, test_tree in enumerate(test_trees):
            expected_tree = expected_trees[idx]
            test_tree.reroot_at_midpoint(update_splits=True)
            self.assertEqual(test_tree.symmetric_difference(expected_tree), 0)
            for split in test_tree.split_edges:
                if test_tree.split_edges[split].head_node is test_tree.seed_node:
                    continue
                self.assertAlmostEqual(test_tree.split_edges[split].length, expected_tree.split_edges[split].length, 3)

if __name__ == "__main__":
    unittest.main()


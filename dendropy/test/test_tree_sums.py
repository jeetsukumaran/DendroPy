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

class TestTreeEdgeSummarization(unittest.TestCase):

    def setUp(self):
        self.taxon_set = dendropy.TaxonSet()
        self.support_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path("primates.beast-mcmc.trees"),
                "nexus",
                taxon_set=self.taxon_set,
                tree_offset=40)
        self.split_distribution = treesplit.SplitDistribution(taxon_set=self.taxon_set)
        self.split_distribution.is_rooted = True
        self.split_distribution.ignore_node_ages = False
        for tree in self.support_trees:
            tree.update_splits()
            self.split_distribution.count_splits_on_tree(tree)

    def testMeanNodeAgeSummarizationOnMCCT(self):
        path_to_target = pathmap.tree_source_path("primates.beast-mcct.noedgelens.tree")
        obs_tree = dendropy.Tree.get_from_path(path_to_target, "nexus")
        obs_tree.update_splits()
        ts = treesum.TreeSummarizer(support_as_labels=True,
                support_as_percentages=False,
                support_label_decimals=4)
        ts.summarize_node_ages_on_tree(tree=obs_tree,
                split_distribution=self.split_distribution,
                set_edge_lengths=True,
                set_extended_attr=True,
                summarization_func=None)
        obs_tree.calc_node_ages()
        exp_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path("primates.beast-mcct.tree"),
                "nexus",
                taxon_set=self.taxon_set)
        exp_tree.update_splits()
        exp_tree.calc_node_ages()
        self.assertEqual(exp_tree.split_edges.keys(), obs_tree.split_edges.keys())
        splits = exp_tree.split_edges.keys()
        for split in splits:
            exp_edge = exp_tree.split_edges[split]
            obs_edge = obs_tree.split_edges[split]
            self.assertAlmostEqual(obs_edge.head_node.age, exp_edge.head_node.age)

if __name__ == "__main__":
    unittest.main()

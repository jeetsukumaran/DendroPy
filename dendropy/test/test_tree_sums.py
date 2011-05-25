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
from dendropy.utility import statistics
from dendropy.test.support import extendedtest

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
        self.support_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path("primates.beast.mcmc.trees"),
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
        path_to_target = pathmap.tree_source_path("primates.beast.mcct.noedgelens.tree")
        obs_tree = dendropy.Tree.get_from_path(path_to_target, "nexus")
        obs_tree.update_splits()
        ts = treesum.TreeSummarizer(support_as_labels=True,
                support_as_percentages=False,
                support_label_decimals=4)
        ts.summarize_node_ages_on_tree(tree=obs_tree,
                split_distribution=self.split_distribution,
                set_edge_lengths=True,
                summarization_func=statistics.median)
        obs_tree.calc_node_ages()
        exp_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path("primates.beast.mcct.medianh.tre"),
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

class TestTopologyCounter(extendedtest.ExtendedTestCase):

    def testSimple(self):
        taxa = dendropy.TaxonSet()
        tree1_str = "[&U] (A,(B,(C,(D,E))));"
        tree2_str = "[&U] (B,(C,(D,(A,E))));"
        tree3_str = "[&U] (D,(A,(B,(C,E))));"
        tree4_str = "[&U] (C,(D,(A,(B,E))));"
        tree5_str = "[&U] (A,(E,(B,(C,D))));"
        all_tree_strs = [tree1_str, tree2_str, tree3_str, tree4_str, tree5_str]
        weights = [8, 5, 4, 2, 1]
        test_tree_strs = []
        for idx, tree_str in enumerate(all_tree_strs):
            test_tree_strs.extend([tree_str] * weights[idx])
        test_trees = dendropy.TreeList.get_from_string("\n".join(test_tree_strs),
                'newick',
                taxon_set=taxa)
        tc = treesum.TopologyCounter()
        expected_freq_values = [float(i)/sum(weights) for i in weights]
        expected_trees = dendropy.TreeList.get_from_string("\n".join(all_tree_strs),
                'newick',
                taxon_set=taxa)
        for tree in test_trees:
            tc.count(tree)
        result_tree_freqs = tc.calc_tree_freqs(taxon_set=taxa)
        for idx, (result_tree, result_freq) in enumerate(result_tree_freqs.items()):
            expected_tree = expected_trees[idx]
            expected_tree.update_splits()
            expected_freq = expected_freq_values[idx]
            expected_count = weights[idx]
            self.assertEqual(result_tree.symmetric_difference(expected_tree), 0,
                    "%s != %s" % (result_tree.as_string('newick'), expected_tree.as_string('newick')))
            self.assertAlmostEqual(result_freq[0], expected_count)
            self.assertAlmostEqual(result_freq[1], expected_freq)

if __name__ == "__main__":
    unittest.main()

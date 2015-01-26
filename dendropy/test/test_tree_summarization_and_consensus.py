#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
from dendropy.calculate import treecompare
from dendropy.test.support import pathmap
from dendropy.mathlib import statistics
from dendropy.test.support import dendropytest

class TestConsensusTree(unittest.TestCase):

    def setUp(self):
        self.tree_list = dendropy.TreeList()
        for t in range(1, 5):
            tf = pathmap.tree_source_path('pythonidae.mb.run%d.t' % t)
            self.tree_list.read_from_path(tf,
                    'nexus',
                    collection_offset=0,
                    tree_offset=25)
        self.mb_con_tree = dendropy.Tree.get_from_path(
                pathmap.tree_source_path("pythonidae.mb.con"),
                schema="nexus",
                taxon_namespace=self.tree_list.taxon_namespace)
        self.mb_con_tree.encode_bipartitions()

    def testConsensus(self):
        con_tree = self.tree_list.consensus(
                min_freq=0.50,
                is_bipartitions_updated=False,
                support_label_decimals=2)
        con_tree.encode_bipartitions()
        self.assertEqual(treecompare.symmetric_difference(self.mb_con_tree, con_tree), 0)
        self.assertEqual(len(con_tree.bipartition_encoding), len(self.mb_con_tree.bipartition_encoding))
        for bipartition in self.mb_con_tree.bipartition_encoding:
            edge1 = self.mb_con_tree.bipartition_edge_map[bipartition]
            edge2 = con_tree.bipartition_edge_map[bipartition]
            if edge1.head_node.label and edge2.head_node.label:
                s1 = float(edge1.head_node.label)
                s2 = round(float(edge2.head_node.label), 2)
                self.assertAlmostEqual(s1, s2, 2)

# class TestTreeCredibilityScoring(unittest.TestCase):

#     def test_from_trees_noburnin_max_product_cc(self):
#         trees = dendropy.TreeList.get_from_path(
#                 pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
#                 "nexus")
#         tsum = treesum.TreeSummarizer()
#         t1, t2 = tsum.calculate_tree_clade_credibilities(trees=trees)
#         for t in trees:
#             self.assertTrue(hasattr(t, "log_product_of_split_support"))
#             self.assertTrue(hasattr(t, "sum_of_split_posteriors"))
#         self.assertEqual(trees.index(t1), 70)
#         self.assertAlmostEqual(t1.log_product_of_split_support, -33.888380488585284)
#         # self.assertAlmostEqual(t1.sum_of_split_posteriors, 85.83000000000001)

#     def test_from_trees_noburnin_max_sum_cc(self):
#         trees = dendropy.TreeList.get_from_path(
#                 pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
#                 "nexus")
#         tsum = treesum.TreeSummarizer()
#         t1, t2 = tsum.calculate_tree_clade_credibilities(trees=trees)
#         for t in trees:
#             self.assertTrue(hasattr(t, "log_product_of_split_support"))
#             self.assertTrue(hasattr(t, "sum_of_split_posteriors"))
#         self.assertEqual(trees.index(t2), 73)
#         # self.assertAlmostEqual(t2.log_product_of_split_support, -38.45253940270466)
#         self.assertAlmostEqual(t2.sum_of_split_posteriors, 89.89000000000001)

#     def test_from_trees_with_burnin_max_product_cc(self):
#         trees = dendropy.TreeList.get_from_path(
#                 pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
#                 "nexus")
#         burnin = 30
#         tsum = treesum.TreeSummarizer()
#         t1, t2 = tsum.calculate_tree_clade_credibilities(
#                 trees=trees,
#                 burnin=burnin)
#         for t in trees[burnin:]:
#             self.assertTrue(hasattr(t, "log_product_of_split_support"))
#             self.assertTrue(hasattr(t, "sum_of_split_posteriors"))

#         # Best tree: bootrep71 (tree number 71)
#         # Highest Log Clade Credibility: -33.95771606695942
#         self.assertEqual(trees.index(t1), 70)
#         self.assertAlmostEqual(t1.log_product_of_split_support, -33.95771606695942)
#         # self.assertAlmostEqual(t1.sum_of_split_posteriors, 85.98571428571427)

#     def test_from_trees_with_burnin_max_sum_cc(self):
#         trees = dendropy.TreeList.get_from_path(
#                 pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
#                 "nexus")
#         burnin = 30
#         tsum = treesum.TreeSummarizer()
#         t1, t2 = tsum.calculate_tree_clade_credibilities(
#                 trees=trees,
#                 burnin=burnin)
#         for t in trees[burnin:]:
#             self.assertTrue(hasattr(t, "log_product_of_split_support"))
#             self.assertTrue(hasattr(t, "sum_of_split_posteriors"))

#         # Best tree: bootrep74 (tree number 74)
#         # Highest Sum Clade Credibility: 30.89
#         # Best tree: bootrep74 (tree number 74)
#         # Highest Sum Clade Credibility: 31.185714285714287
#         self.assertTrue(trees.index(t2), 73)
#         # self.assertAlmostEqual(t2.log_product_of_split_support, -37.912350577390605)
#         self.assertAlmostEqual(t2.sum_of_split_posteriors, 31.185714285714287)

class TestTreeEdgeSummarization(unittest.TestCase):

    def setUp(self):
        self.support_trees_path = pathmap.tree_source_path("primates.beast.mcmc.trees")
        self.target_tree_path = pathmap.tree_source_path("primates.beast.mcct.noedgelens.tree")
        self.expected_tree_path = pathmap.tree_source_path("primates.beast.mcct.medianh.tre")
        self.burnin = 40

    def testMeanNodeAgeSummarizationOnMCCT(self):
        tree_array = dendropy.TreeArray(ignore_node_ages=False)
        tree_array.read_from_path(
                self.support_trees_path,
                "nexus",
                # colleciton_offset=0,
                tree_offset=self.burnin,
                )
        target_tree = dendropy.Tree.get_from_path(
                self.target_tree_path,
                schema="nexus",
                taxon_namespace=tree_array.taxon_namespace,
                )
        tree_array.summarize_splits_on_tree(
                tree=target_tree,
                is_bipartitions_updated=False,
                set_edge_lengths="median-age",
                )
        expected_tree = dendropy.Tree.get_from_path(
                self.expected_tree_path,
                "nexus",
                taxon_namespace=tree_array.taxon_namespace)
        expected_tree.encode_bipartitions()
        expected_tree.calc_node_ages()
        self.assertEqual(expected_tree.bipartition_encoding, target_tree.bipartition_encoding)
        for exp_bipartition in expected_tree.bipartition_encoding:
            exp_edge = expected_tree.bipartition_edge_map[exp_bipartition]
            obs_edge = target_tree.bipartition_edge_map[exp_bipartition]
            self.assertAlmostEqual(obs_edge.head_node.age, exp_edge.head_node.age)

class TestTopologyCounter(dendropytest.ExtendedTestCase):

    def get_regime(self, is_rooted, weights):
        tree_strings = [
            "(A,(B,(C,(D,E))));",
            "(B,(C,(D,(A,E))));",
            "(D,(A,(B,(C,E))));",
            "(C,(D,(A,(B,E))));",
            "(A,(E,(B,(C,D))));",
            ]
        if is_rooted is None:
            rooting = ""
        if is_rooted:
            rooting = " [&R] "
        else:
            rooting = " [&U] "
        if weights is None:
            weightings = ["" for w in len(tree_strings)]
        else:
            assert len(weights) == len(tree_strings)
            weightings = [" [&W {}] ".format(weight) for w in weights]
        for idx, (tree_string, weight) in zip(tree_strings, weightings):
            tree_strings[idx] = "{}{}{}".format(rooting, weight, tree_string)
        self.taxon_namespace = dendropy.TaxonNamespace()

    def testSimple(self):
        self.taxon_namespace = dendropy.TaxonNamespace()
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
        test_trees = dendropy.TreeList.get_from_string(
                "\n".join(test_tree_strs),
                'newick',
                taxon_namespace=self.taxon_namespace)
        # expected_freq_values = [float(i)/sum(weights) for i in weights]
        expected_trees = dendropy.TreeList.get_from_string(
                "\n".join(all_tree_strs),
                'newick',
                taxon_namespace=self.taxon_namespace)
        expected_freqs = {}
        for idx, tree in enumerate(expected_trees):
            b = frozenset(tree.encode_bipartitions())
            expected_freqs[b] = float(weights[idx])/sum(weights)
        ta = test_trees.as_tree_array()
        result_tree_freqs = ta.topology_frequencies()
        self.assertEqual(len(result_tree_freqs), len(expected_freqs))
        for tree in result_tree_freqs:
            b = frozenset(tree.encode_bipartitions())
            self.assertAlmostEqual(result_tree_freqs[tree], expected_freqs[b])

if __name__ == "__main__":
    unittest.main()

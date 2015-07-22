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
Tests of summarization.
"""

import collections
import unittest
import dendropy
import random
import itertools
from dendropy.calculate import treecompare
from dendropy.test.support import pathmap
from dendropy.calculate import statistics
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

class TestBasicCredibilityScoring(unittest.TestCase):

    def get_trees(self):
        trees = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
                "nexus")
        return trees

    def setUp(self):
        self.trees = self.get_trees()

    def test_product_of_credibilities(self):
        ta = self.trees.as_tree_array(is_rooted_trees=True)
        sd = self.get_trees().split_distribution(is_bipartitions_updated=False) # for independent verification
        scores, max_idx = ta.calculate_log_product_of_split_supports()
        self.assertEqual(len(scores), len(self.trees))
        for score, tree in zip(scores, self.trees):
            self.assertAlmostEqual(score, sd.log_product_of_split_support_on_tree(tree))
        self.assertEqual(max_idx, 70)
        self.assertAlmostEqual(scores[max_idx], -33.888380488585284)
        t0 = self.trees[70]
        t1 = ta.maximum_product_of_split_support_tree()
        self.assertEqual(treecompare.symmetric_difference(t0, t1), 0)

    def test_sum_of_credibilities(self):
        ta = self.trees.as_tree_array(is_rooted_trees=True)
        sd = self.get_trees().split_distribution(is_bipartitions_updated=False) # for independent verification
        scores, max_idx = ta.calculate_sum_of_split_supports()
        self.assertEqual(len(scores), len(self.trees))
        for score, tree in zip(scores, self.trees):
            self.assertAlmostEqual(score, sd.sum_of_split_support_on_tree(tree))
        self.assertEqual(max_idx, 73)
        self.assertAlmostEqual(scores[max_idx], 30.89)
        t0 = self.trees[73]
        t1 = ta.maximum_sum_of_split_support_tree()
        self.assertEqual(treecompare.symmetric_difference(t0, t1), 0)

    def test_split_distribution_max_sum_of_credibilities(self):
        sd = self.trees.split_distribution(is_bipartitions_updated=False)
        t0 = self.trees[73]
        score = sd.sum_of_split_support_on_tree(t0)
        self.assertAlmostEqual(score, 30.89)

        # # Best tree: bootrep74 (tree number 74)
        # # Highest Sum Clade Credibility: 30.89
        # # Best tree: bootrep74 (tree number 74)
        # # Highest Sum Clade Credibility: 31.185714285714287

        # scores, max_idx = ta.calculate_log_product_of_split_supports()
        # self.assertEqual(len(scores), len(self.trees))
        # self.assertEqual(max_idx, 70)
        # t1 = ta.maximum_product_of_split_support_tree()
        # self.assertEqual(treecompare.symmetric_difference(t0, t1), 0)
        # t1, t2 = tsum.calculate_tree_clade_credibilities(trees=trees)
        # for t in trees:
        #     self.assertTrue(hasattr(t, "log_product_of_split_support"))
        #     self.assertTrue(hasattr(t, "sum_of_split_posteriors"))
        # self.assertEqual(trees.index(t1), 70)
        # self.assertAlmostEqual(t1.log_product_of_split_support, -33.888380488585284)
        # # self.assertAlmostEqual(t1.sum_of_split_posteriors, 85.83000000000001)

    # def test_from_trees_noburnin_max_sum_cc(self):
        # trees = dendropy.TreeList.get_from_path(
        #         pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
        #         "nexus")
        # tsum = treesum.TreeSummarizer()
        # t1, t2 = tsum.calculate_tree_clade_credibilities(trees=trees)
        # for t in trees:
        #     self.assertTrue(hasattr(t, "log_product_of_split_support"))
        #     self.assertTrue(hasattr(t, "sum_of_split_posteriors"))
        # self.assertEqual(trees.index(t2), 73)
        # # self.assertAlmostEqual(t2.log_product_of_split_support, -38.45253940270466)
        # self.assertAlmostEqual(t2.sum_of_split_posteriors, 89.89000000000001)

    # def test_from_trees_with_burnin_max_product_cc(self):
        # trees = dendropy.TreeList.get_from_path(
        #         pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
        #         "nexus")
        # burnin = 30
        # tsum = treesum.TreeSummarizer()
        # t1, t2 = tsum.calculate_tree_clade_credibilities(
        #         trees=trees,
        #         burnin=burnin)
        # for t in trees[burnin:]:
        #     self.assertTrue(hasattr(t, "log_product_of_split_support"))
        #     self.assertTrue(hasattr(t, "sum_of_split_posteriors"))

        # # Best tree: bootrep71 (tree number 71)
        # # Highest Log Clade Credibility: -33.95771606695942
        # self.assertEqual(trees.index(t1), 70)
        # self.assertAlmostEqual(t1.log_product_of_split_support, -33.95771606695942)
        # # self.assertAlmostEqual(t1.sum_of_split_posteriors, 85.98571428571427)

    # def test_from_trees_with_burnin_max_sum_cc(self):
        # trees = dendropy.TreeList.get_from_path(
        #         pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
        #         "nexus")
        # burnin = 30
        # tsum = treesum.TreeSummarizer()
        # t1, t2 = tsum.calculate_tree_clade_credibilities(
        #         trees=trees,
        #         burnin=burnin)
        # for t in trees[burnin:]:
        #     self.assertTrue(hasattr(t, "log_product_of_split_support"))
        #     self.assertTrue(hasattr(t, "sum_of_split_posteriors"))

        # # Best tree: bootrep74 (tree number 74)
        # # Highest Sum Clade Credibility: 30.89
        # # Best tree: bootrep74 (tree number 74)
        # # Highest Sum Clade Credibility: 31.185714285714287
        # self.assertTrue(trees.index(t2), 73)
        # # self.assertAlmostEqual(t2.log_product_of_split_support, -37.912350577390605)
        # self.assertAlmostEqual(t2.sum_of_split_posteriors, 31.185714285714287)

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

    def get_regime(self,
            is_rooted,
            is_multifurcating,
            is_weighted,
            tree_offset=0,
            taxon_namespace=None,
            num_trees=500):
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        if is_multifurcating:
            if is_rooted:
                tree_filename = "dendropy-test-trees-multifurcating-rooted.nexus"
            else:
                tree_filename = "dendropy-test-trees-multifurcating-unrooted.nexus"
        else:
            if is_rooted:
                tree_filename = "dendropy-test-trees-n10-rooted-treeshapes.nexus"
            else:
                tree_filename = "dendropy-test-trees-n14-unrooted-treeshapes.nexus"
        source_trees = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path(tree_filename),
                "nexus",
                taxon_namespace=taxon_namespace)
        for tree in source_trees:
            tree.encode_bipartitions()
            tree.key = frozenset(tree.bipartition_encoding)
            tree.total_weighted_count = 0.0
            tree.actual_count = 0
        # if is_weighted:
        #     weights = []
        #     for tree in source_trees:
        #         w = random.uniform(0.1, 10)
        #         tree.weight = w
        #         weights.append(w)
        # else:
        #     weights = [1.0 for i in len(source_trees)]
        test_tree_strings = []
        total_weight = 0.0
        while len(test_tree_strings) < num_trees:
            tree = random.choice(source_trees)
            if len(test_tree_strings) >= tree_offset:
                tree.actual_count += 1
            if is_weighted:
                weight = random.choice([0.25, 1.0, 2.8, 5.6, 11.0,])
                tree.weight = weight
                if len(test_tree_strings) >= tree_offset:
                    tree.total_weighted_count += weight
                    total_weight += weight
            else:
                tree.weight = None
                if len(test_tree_strings) >= tree_offset:
                    tree.total_weighted_count += 1.0
                    total_weight += 1.0
            for nd in tree:
                nd.edge.length = random.uniform(0, 100)
            test_tree_strings.append(tree.as_string(
                schema="newick",
                store_tree_weights=is_weighted,
                suppress_edge_lengths=False,
                suppress_internal_node_labels=True,
                suppress_internal_taxon_labels=True,
                ))
        test_trees_string = "\n".join(test_tree_strings)
        bipartition_encoding_freqs = {}
        source_trees.total_weight = total_weight
        for tree in source_trees:
            tree.frequency = float(tree.total_weighted_count) / total_weight
            bipartition_encoding_freqs[tree.key] = tree.frequency
        return source_trees, bipartition_encoding_freqs, test_trees_string

    def testVariants(self):
        for tree_offset, is_weighted, is_multifurcating, is_rooted in itertools.product( (100,), (False, True, ), (False, True, ), (False, True, ),  ):
        # for tree_offset, is_weighted, is_multifurcating, is_rooted in itertools.product( (0, 100), (True,), (False,), (False,),  ):
            # print("is_rooted: {is_rooted}, is_multifurcating: {is_multifurcating}, is_weighted: {is_weighted}, tree_offset: {tree_offset}".format(
            #     is_rooted=is_rooted,
            #     is_multifurcating=is_multifurcating,
            #     is_weighted=is_weighted,
            #     tree_offset=tree_offset))
            source_trees, bipartition_encoding_freqs, test_trees_string = self.get_regime(
                    is_rooted=is_rooted,
                    is_multifurcating=is_multifurcating,
                    is_weighted=is_weighted,
                    tree_offset=tree_offset)
            ta = dendropy.TreeArray(
                    is_rooted_trees=is_rooted,
                    use_tree_weights=is_weighted,
                    taxon_namespace=source_trees.taxon_namespace,
                    )
            ta.read_from_string(
                    test_trees_string,
                    "newick",
                    tree_offset=tree_offset,
                    store_tree_weights=is_weighted)
            be_to_tree = {}
            for tree in source_trees:
                be_to_tree[tree.key] = tree
            topologies = ta.topologies()
            for tree in topologies:
                b = frozenset(tree.encode_bipartitions())
                # stree = be_to_tree[b]
                # print("{} ({}): {}".format(
                #     calculated_topology_freqs[tree],
                #     ta._split_distribution.calc_normalization_weight(),
                #     (   bipartition_encoding_freqs[b],
                #         stree.actual_count,
                #         stree.total_weighted_count,
                #         source_trees.total_weight,
                #         stree.frequency,
                #         stree.total_weighted_count / source_trees.total_weight,
                #     )))
                self.assertAlmostEqual(
                        tree.frequency,
                        bipartition_encoding_freqs[b])

            calculated_bipartition_encoding_freqs = ta.bipartition_encoding_frequencies()
            for tree in source_trees:
                # if tree.key not in calculated_bipartition_encoding_freqs:
                #     print(tree.actual_count)
                #     print(tree.total_weighted_count)
                #     print(tree.frequency)
                # f1 = bipartition_encoding_freqs[tree.key]
                # f2 = calculated_bipartition_encoding_freqs[tree.key]
                # self.assertAlmostEqual(f1,f2)
                if tree.actual_count == 0:
                    if tree.key in calculated_bipartition_encoding_freqs:
                        self.assertAlmostEqual(calculated_bipartition_encoding_freqs[tree.key], 0)
                else:
                    # self.assertIn(tree.key, calculated_bipartition_encoding_freqs)
                    f1 = bipartition_encoding_freqs[tree.key]
                    f2 = calculated_bipartition_encoding_freqs[tree.key]
                    self.assertAlmostEqual(f1,f2)

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
        topologies = ta.topologies()
        self.assertEqual(len(topologies), len(expected_freqs))
        for tree in topologies:
            b = frozenset(tree.encode_bipartitions())
            self.assertAlmostEqual(tree.frequency, expected_freqs[b])

if __name__ == "__main__":
    unittest.main()

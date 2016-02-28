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
Testing of calculation of and operations with split bitmask hashes.
"""

import warnings
import unittest
import re
import sys
from dendropy.test.support import pathmap
from dendropy.test.support import paupsplitsreference
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.utility import messaging
from dendropy.utility import bitprocessing
from dendropy.interop import paup
from dendropy.utility.textprocessing import StringIO
from dendropy.calculate import treecompare
import dendropy

_LOG = messaging.get_logger(__name__)

class SplitDistributionTestCases(ExtendedTestCase):

    def check_splits_distribution(self,
            tree_filename,
            splits_filename,
            use_tree_weights,
            is_rooted,
            expected_num_trees,
            ):
        if is_rooted is None:
            key_column_index = 2 # default to unrooted: normalized split bitmask
        elif is_rooted:
            key_column_index = 1 # leafset_bitmask / unnormalized split bitmask
        else:
            key_column_index = 2 # normalized split bitmask
        splits_ref = paupsplitsreference.get_splits_reference(
                splits_filename=splits_filename,
                key_column_index=key_column_index,
                )
        # print("* {} ({})".format(tree_filename, splits_filename))
        tree_filepath = pathmap.tree_source_path(tree_filename)
        trees = dendropy.TreeList.get_from_path(
                tree_filepath,
                "nexus",
                store_tree_weights=use_tree_weights)
        sd = dendropy.SplitDistribution(
                taxon_namespace=trees.taxon_namespace,
                use_tree_weights=use_tree_weights)
        for tree in trees:
            sd.count_splits_on_tree(tree)

        # trees counted ...
        self.assertEqual(sd.total_trees_counted, len(trees))
        # frequencies have not yet been calculated
        self.assertEqual(sd._trees_counted_for_freqs, 0)
        self.assertFalse(sd.is_mixed_rootings_counted())
        if is_rooted:
            self.assertTrue(sd.is_all_counted_trees_rooted())
        else:
            self.assertFalse(sd.is_all_counted_trees_rooted())
            self.assertTrue(sd.is_all_counted_trees_treated_as_unrooted() or sd.is_all_counted_trees_strictly_unrooted())

        # splits_distribution also counts trivial splits, so this will not work
        # self.assertEqual(len(splits_ref), len(sd))

        expected_nontrivial_splits = list(splits_ref.keys())
        observed_splits = set(sd.split_counts.keys())
        visited_splits = []
        # for k in sorted(observed_splits):
        #     print("{}: {}, {}".format(k, sd.split_counts[k], sd[k]))
        all_taxa_bitmask = sd.taxon_namespace.all_taxa_bitmask()
        for split in expected_nontrivial_splits:
            self.assertAlmostEqual(sd.split_counts[split], splits_ref[split]["count"], 2,
                    "{} (using '{}'): {}".format(tree_filename, splits_filename, split))
            self.assertAlmostEqual(sd[split], splits_ref[split]["frequency"], 2,
                    "{} (using '{}'): {}".format(tree_filename, splits_filename, split))
            self.assertAlmostEqual(sd.split_frequencies[split], splits_ref[split]["frequency"], 2,
                    "{} (using '{}'): {}".format(tree_filename, splits_filename, split))
            if split in observed_splits:
                observed_splits.remove(split)
            visited_splits.append(split)
        self.assertEqual(len(visited_splits), len(expected_nontrivial_splits))

        # ensure remaining splits (not given in PAUP splits file) are trivial ones (which are not tracked by PAUP)
        for split in observed_splits:
            self.assertTrue(dendropy.Bipartition.is_trivial_bitmask(split, all_taxa_bitmask))

    def test_group1(self):
        sources = [
                ("cetaceans.mb.no-clock.mcmc.trees"    , 251, False, False), # Trees explicitly unrooted
                ("cetaceans.mb.no-clock.mcmc.weighted-01.trees" , 251, False , True), # Weighted
                ("cetaceans.mb.no-clock.mcmc.weighted-02.trees" , 251, False , True), # Weighted
                ("cetaceans.mb.no-clock.mcmc.weighted-03.trees" , 251, False , True), # Weighted
                ("cetaceans.mb.strict-clock.mcmc.trees", 251, True , False), # Trees explicitly rooted
                ("cetaceans.mb.strict-clock.mcmc.weighted-01.trees" , 251, True , True), # Weighted
                ("cetaceans.mb.strict-clock.mcmc.weighted-02.trees" , 251, True , True), # Weighted
                ("cetaceans.mb.strict-clock.mcmc.weighted-03.trees" , 251, True , True), # Weighted
                ("issue_mth_2009-02-03.rooted.nexus"   , 100, True , False), # 100 trees (frequency column not reported by PAUP)
                ("issue_mth_2009-02-03.unrooted.nexus" , 100, False , False), # 100 trees (frequency column not reported by PAUP)
                ("cetaceans.raxml.bootstraps.trees"    , 250, None , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                ("cetaceans.raxml.bootstraps.weighted-01.trees"    , 250, None , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                ("cetaceans.raxml.bootstraps.weighted-02.trees"    , 250, None , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                ("cetaceans.raxml.bootstraps.weighted-03.trees"    , 250, None , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
        ]
        splits_filename_template = "{stemname}.is-rooted-{is_rooted}.use-tree-weights-{use_weights}.burnin-{burnin}.splits.txt"
        for tree_filename, num_trees, treefile_is_rooted, treefile_is_weighted in sources:
            stemname = tree_filename
            for use_weights in (False, True, None):
                expected_is_rooted = treefile_is_rooted
                splits_filename = splits_filename_template.format(
                        stemname=stemname,
                        is_rooted=expected_is_rooted,
                        use_weights=use_weights,
                        burnin=0)
                self.check_splits_distribution(
                        tree_filename=tree_filename,
                        splits_filename=splits_filename,
                        is_rooted=treefile_is_rooted,
                        use_tree_weights=use_weights,
                        expected_num_trees=num_trees)

if not paup.DENDROPY_PAUP_INTEROPERABILITY:
    _LOG.warn("PAUP interoperability not available: skipping split counting tests")
else:

    class SplitCountTest(ExtendedTestCase):

        @classmethod
        def setUpClass(cls):
            if sys.version_info.major < 3:
                cls.assertRaisesRegex = cls.assertRaisesRegexp

        def check_split_counting(self,
                tree_filename,
                test_as_rooted,
                parser_rooting_interpretation,
                test_ignore_tree_weights=False,
                dp_ignore_tree_weights=False,
                ):
            tree_filepath = pathmap.tree_source_path(tree_filename)
            ps = paup.PaupService()
            paup_sd = ps.get_split_distribution_from_files(
                    tree_filepaths=[tree_filepath],
                    is_rooted=test_as_rooted,
                    use_tree_weights=not test_ignore_tree_weights,
                    burnin=0,
                    taxa_definition_filepath=tree_filepath
                    )
            taxon_namespace = paup_sd.taxon_namespace
            dp_sd = dendropy.SplitDistribution(taxon_namespace=taxon_namespace)
            dp_sd.ignore_edge_lengths = True
            dp_sd.ignore_node_ages = True
            dp_sd.ignore_tree_weights = dp_ignore_tree_weights
            taxa_mask = taxon_namespace.all_taxa_bitmask()
            taxon_namespace.is_mutable = False
            trees = dendropy.TreeList.get_from_path(tree_filepath,
                    "nexus",
                    rooting=parser_rooting_interpretation,
                    taxon_namespace=taxon_namespace)
            for tree in trees:
                self.assertIs(tree.taxon_namespace, taxon_namespace)
                self.assertIs(tree.taxon_namespace, dp_sd.taxon_namespace)
                dp_sd.count_splits_on_tree(
                        tree,
                        is_bipartitions_updated=False)
            self.assertEqual(dp_sd.total_trees_counted, paup_sd.total_trees_counted)
            taxa_mask = taxon_namespace.all_taxa_bitmask()
            for split in dp_sd.split_counts:
                if not dendropy.Bipartition.is_trivial_bitmask(split, taxa_mask):
                    # if split not in paup_sd.split_counts:
                    #     print("{}: {}".format(split, split in paup_sd.split_counts))
                    #     s2 = taxon_namespace.normalize_bitmask(split)
                    #     print("{}: {}".format(s2, s2 in paup_sd.split_counts))
                    #     s3 = ~split & taxon_namespace.all_taxa_bitmask()
                    #     print("{}: {}".format(s3, s3 in paup_sd.split_counts))
                    self.assertIn(split, paup_sd.split_counts, "split not found")
                    self.assertEqual(dp_sd.split_counts[split], paup_sd.split_counts[split], "incorrect split frequency")
                    del paup_sd.split_counts[split]
            remaining_splits = list(paup_sd.split_counts.keys())
            for split in remaining_splits:
                if dendropy.Bipartition.is_trivial_bitmask(split, taxa_mask):
                    del paup_sd.split_counts[split]
            self.assertEqual(len(paup_sd.split_counts), 0)

        def test_basic_split_count_with_incorrect_rootings_raises_error(self):
            assertion_error_regexp1 = re.compile("(incorrect split frequency|split not found)")
            test_cases = (
                ('pythonidae.reference-trees.nexus', True, "force-unrooted", assertion_error_regexp1),
                ('feb032009.trees.nexus', False, "force-rooted", assertion_error_regexp1),
                )
            for test_case, test_as_rooted, parser_rooting_interpretation, assertion_error_regexp in test_cases:
                with self.assertRaisesRegex(AssertionError, assertion_error_regexp):
                    self.check_split_counting(
                            test_case,
                            test_as_rooted=test_as_rooted,
                            parser_rooting_interpretation=parser_rooting_interpretation)

        def test_basic_split_count_with_incorrect_weight_treatment_raises_error(self):
            assertion_error_regexp1 = re.compile("incorrect split frequency")
            test_cases = (
                    ("cetaceans.mb.no-clock.mcmc.weighted-01.trees", False),
                    ("cetaceans.mb.strict-clock.mcmc.weighted-01.trees", True),
                )
            for test_case, test_as_rooted in test_cases:
                with self.assertRaisesRegex(AssertionError, assertion_error_regexp1):
                    self.check_split_counting(
                            test_case,
                            test_as_rooted=test_as_rooted,
                            parser_rooting_interpretation="default-rooted",
                            test_ignore_tree_weights=False,
                            dp_ignore_tree_weights=False,
                            )

        def test_basic_split_counting_under_different_rootings(self):
            test_cases = (
                'pythonidae.reference-trees.nexus',
                'feb032009.trees.nexus',
                'maj-rule-bug1.trees.nexus',
                'maj-rule-bug2.trees.nexus',
                )
            for is_rooted in (True, False):
                if is_rooted:
                    rooting = "force-rooted"
                else:
                    rooting = "force-unrooted"
                for test_case in test_cases:
                    self.check_split_counting(
                            test_case,
                            test_as_rooted=is_rooted,
                            parser_rooting_interpretation=rooting)

class CladeMaskTest(unittest.TestCase):

    def runTest(self):
        # rooted tree: so clade bitmasks
        tree_list = dendropy.TreeList.get_from_stream(
            StringIO("""[&R]((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);"""),
            "newick")
        for i in tree_list:
            _LOG.debug(i._get_indented_form())
            i.encode_bipartitions()
            _LOG.debug(i._get_indented_form(splits=True))
            i._debug_check_tree(splits=True, logger_obj=_LOG)
        root1 = tree_list[0].seed_node
        root1e = root1.edge
        self.assertEqual(bitprocessing.indexes_of_set_bits(root1e.split_bitmask), list(range(6)))
        self.assertEqual(bitprocessing.indexes_of_set_bits(root1e.split_bitmask, one_based=True), list(range(1,7)))
        self.assertEqual(bitprocessing.indexes_of_set_bits(root1e.split_bitmask, fill_bitmask=21, one_based=True), [1, 3, 5])
        self.assertEqual(bitprocessing.indexes_of_set_bits(root1e.split_bitmask, fill_bitmask=21), [0, 2, 4])
        self.assertEqual(bitprocessing.num_set_bits(root1e.split_bitmask), 6)

        fc1 = root1.child_nodes()[0]
        fc1e = fc1.edge
        self.assertEqual(bitprocessing.indexes_of_set_bits(fc1e.split_bitmask), [0, 1])
        self.assertEqual(bitprocessing.indexes_of_set_bits(fc1e.split_bitmask, one_based=True), [1, 2])
        self.assertEqual(bitprocessing.indexes_of_set_bits(fc1e.split_bitmask, fill_bitmask=0x15, one_based=True), [1])
        self.assertEqual(bitprocessing.indexes_of_set_bits(fc1e.split_bitmask, fill_bitmask=0x15), [0])
        self.assertEqual(bitprocessing.num_set_bits(fc1e.split_bitmask), 2)

class CountBitsTest(unittest.TestCase):

    def runTest(self):
        self.assertEqual(bitprocessing.num_set_bits(21), 3)

class LowestBitTest(unittest.TestCase):

    def runTest(self):
        for n, expected in enumerate([0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16]):
            self.assertEqual(bitprocessing.least_significant_set_bit(n), expected)

class IsTrivialTest(unittest.TestCase):

    def runTest(self):
        y = True
        n = False
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, dendropy.Bipartition.is_trivial_bitmask(i, 0xF))
        for i, r in enumerate([y, y, y, n, y, n, n, n, y, n, n, n, n, n, n, y, y, n, n, n, n, n, n, y, n, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, dendropy.Bipartition.is_trivial_bitmask(i, 0x1F))
                              #0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, dendropy.Bipartition.is_trivial_bitmask(i, 0x17))

class IncompleteLeafSetSplitTest(unittest.TestCase):

    def check(self, title, src_prefix):
        tns = dendropy.TaxonNamespace()
        input_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".dendropy-pruned.nex"),
                schema='nexus',
                attached_taxon_namespace=tns)
        input_taxa = input_ds.taxon_namespaces[0]
        output_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".paup-pruned.nex"),
                schema='nexus',
                taxon_namespace=input_taxa)
        for set_idx, src_trees in enumerate(input_ds.tree_lists):
            src_trees = input_ds.tree_lists[set_idx]
            ref_trees = output_ds.tree_lists[set_idx]
            for tree_idx, src_tree in enumerate(src_trees):
                _LOG.debug("%s Set %d/%d, Tree %d/%d" % (title, set_idx+1, len(input_ds.tree_lists), tree_idx+1, len(src_trees)))
                ref_tree = ref_trees[tree_idx]
                # tree_dist = paup.symmetric_difference(src_tree, ref_tree)
                # d = src_tree.symmetric_difference(ref_tree)
                # if d > 0:
                #     print d
                self.assertEqual(treecompare.symmetric_difference(src_tree, ref_tree), 0)

    def testUnrooted(self):
        self.check("Unrooted", "incomplete_leaves_unrooted")

    def testRooted(self):
        self.check("Rooted", "incomplete_leaves_rooted")

    def testPrunedThenEncoding(self):
        inp = StringIO('''(a,b,c,(d,e));
        (b,d,(c,e));''')
        first, second = dendropy.TreeList.get_from_stream(inp, schema='newick')
        # prune tree 1 to have the same leaf set as tree 2.
        #   this removes the first taxon in the taxon list "A"
        retain_list = set([node.taxon for node in second.leaf_nodes()])
        exclude_list = [node for node in first.leaf_nodes() if node.taxon not in retain_list]
        for nd in exclude_list:
            first.prune_subtree(nd)
        # the trees are now (b,c,(d,e)) and (b,d,(c,e)) so the symmetric diff is 2
        self.assertEqual(2, treecompare.symmetric_difference(first, second))

class TestTreeSplitSupportCredibilityScoring(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path("issue_mth_2009-02-03.rooted.nexus"),
                "nexus")
        self.split_distribution = dendropy.SplitDistribution(taxon_namespace=self.trees.taxon_namespace)
        for tree in self.trees:
            self.split_distribution.count_splits_on_tree(
                    tree,
                    is_bipartitions_updated=False)

    def test_product_of_split_support_on_tree(self):
        t1 = self.trees[70]
        self.assertAlmostEqual(
                self.split_distribution.log_product_of_split_support_on_tree(t1),
                -33.888380488585284)

    def test_sum_of_split_support_on_tree(self):
        t1 = self.trees[73]
        self.assertAlmostEqual(
                self.split_distribution.sum_of_split_support_on_tree(t1),
                30.89000000000001)

    def test_sum_of_split_support_on_tree2(self):
        t1 = self.trees[73]
        self.assertAlmostEqual(
                self.split_distribution.sum_of_split_support_on_tree(t1, include_external_splits=True),
                30.89000000000001 + len(self.trees.taxon_namespace))

if __name__ == "__main__":
    unittest.main()
    # if paup.DENDROPY_PAUP_INTEROPERABILITY:
    #     unittest.main()
    # else:
    #     _LOG.warn("PAUP interoperability not available: skipping split counting tests")

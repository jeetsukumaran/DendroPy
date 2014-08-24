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
Testing of calculation of and operations with split bitmask hashes.
"""

import warnings
import unittest
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3

from dendropy.test.support import pathmap
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.utility import messaging
from dendropy.interop import paup
from dendropy import treesplit
import dendropy

_LOG = messaging.get_logger(__name__)

class SplitCountTest(ExtendedTestCase):

    def setUp(self):
        self.test_cases = [('pythonidae.reference-trees.nexus', 'pythonidae.reference-trees.nexus')]
        if True: # runlevel.is_test_enabled(runlevel.SLOW, _LOG, self.__class__.__name__):
            self.test_cases.extend([
                ('feb032009.trees.nexus', 'feb032009.trees.nexus'),
                ('maj-rule-bug1.trees.nexus', 'maj-rule-bug1.trees.nexus'),
                ('maj-rule-bug2.trees.nexus', 'maj-rule-bug2.trees.nexus'),
            ])

    def countSplits(self, tc, is_rooted):
        # _LOG.info(tc[0] + "; " + tc[1])
        tree_filepaths = [pathmap.tree_source_path(tc[0])]
        taxa_filepath = pathmap.tree_source_path(tc[1])
        paup_sd = paup.get_split_distribution(
                tree_filepaths,
                taxa_filepath,
                is_rooted=is_rooted,
                burnin=0)
        taxon_namespace = paup_sd.taxon_namespace
        dp_sd = treesplit.SplitDistribution(taxon_namespace=taxon_namespace)
        dp_sd.ignore_edge_lengths = True
        dp_sd.ignore_node_ages = True
        dp_sd.is_rooted = is_rooted

        _LOG.debug("Taxon set: %s" % [t.label for t in taxon_namespace])
        taxa_mask = taxon_namespace.all_taxa_bitmask()
        taxon_namespace.is_mutable = False
        if is_rooted:
            rooting = "force-rooted"
        else:
            rooting = "force-unrooted"
        for tree_filepath in tree_filepaths:
            trees = dendropy.TreeList.get_from_path(tree_filepath,
                    "nexus",
                    rooting=rooting,
                    taxon_namespace=taxon_namespace)
            for tree in trees:
                self.assertIs(tree.taxon_namespace, taxon_namespace)
                self.assertIs(tree.taxon_namespace, dp_sd.taxon_namespace)
                treesplit.encode_splits(tree)
                dp_sd.count_splits_on_tree(tree)
        self.assertEqual(dp_sd.total_trees_counted, paup_sd.total_trees_counted)

        # SplitsDistribution counts trivial splits, whereas PAUP*
        # contree does not, so the following will not work
#            assert len(dp_sd.splits) == len(paup_sd.splits),\
#                 "dp = %d, sd = %d" % (len(dp_sd.splits), len(paup_sd.splits))

        taxa_mask = taxon_namespace.all_taxa_bitmask()
        for split in dp_sd.splits:
            if not treesplit.is_trivial_split(split, taxa_mask):
                self.assertIn(split, paup_sd.splits)
                self.assertEqual(dp_sd.split_counts[split], paup_sd.split_counts[split])
                paup_sd.splits.remove(split)

        # if any splits remain, they were not
        # in dp_sd or were trivial
        remaining_splits = list(paup_sd.splits)
        for split in remaining_splits:
            if treesplit.is_trivial_split(split, taxa_mask):
                paup_sd.splits.remove(split)
        self.assertEqual(len(paup_sd.splits), 0)

    def testUnrootedSplitCounts(self):
        for tc in self.test_cases:
            self.countSplits(tc, is_rooted=False)

    def testRootedSplitCounts(self):
        for tc in self.test_cases:
            self.countSplits(tc, is_rooted=True)

class CladeMaskTest(unittest.TestCase):

    def runTest(self):
        tree_list = dendropy.TreeList.get_from_stream(
            StringIO("""((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);"""),
            "newick")
        for i in tree_list:
            _LOG.debug(i._get_indented_form())
            treesplit.encode_splits(i)
            _LOG.debug(i._get_indented_form(splits=True))
            i._debug_check_tree(splits=True, logger_obj=_LOG)
        root1 = tree_list[0].seed_node
        root1e = root1.edge
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask), list(range(6)))
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, one_based=True), list(range(1,7)))
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, mask=21, one_based=True), [1, 3, 5])
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, mask=21), [0, 2, 4])
        self.assertEqual(treesplit.count_bits(root1e.split_bitmask), 6)

        fc1 = root1.child_nodes()[0]
        fc1e = fc1.edge
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask), [0, 1])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, one_based=True), [1, 2])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, mask=0x15, one_based=True), [1])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, mask=0x15), [0])
        self.assertEqual(treesplit.count_bits(fc1e.split_bitmask), 2)

class CountBitsTest(unittest.TestCase):

    def runTest(self):
        self.assertEqual(treesplit.count_bits(21), 3)

class LowestBitTest(unittest.TestCase):

    def runTest(self):
        for n, expected in enumerate([0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16]):
            self.assertEqual(treesplit.lowest_bit_only(n), expected)

class IsTrivialTest(unittest.TestCase):

    def runTest(self):
        y = True
        n = False
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0xF))
        for i, r in enumerate([y, y, y, n, y, n, n, n, y, n, n, n, n, n, n, y, y, n, n, n, n, n, n, y, n, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0x1F))
                              #0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0x17))

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
                taxon_set=input_taxa)
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
                self.assertEqual(src_tree.symmetric_difference(ref_tree), 0)

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
        self.assertEqual(2, first.symmetric_difference(second))

if __name__ == "__main__":
    if paup.DENDROPY_PAUP_INTEROPERABILITY:
        unittest.main()
    else:
        _LOG.warn("PAUP interoperability not available: skipping split counting tests")

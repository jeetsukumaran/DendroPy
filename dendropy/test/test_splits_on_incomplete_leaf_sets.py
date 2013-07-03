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
Tests of split counting on trees with incomplete leaf sets.
"""

import unittest

import dendropy
from dendropy.test.support import pathmap
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class IncompleteLeafSetSplitTest(unittest.TestCase):

    def check(self, title, src_prefix):
        input_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".dendropy-pruned.nex"),
                schema='nexus',
                attach_taxon_set=True)
        input_taxa = input_ds.taxon_sets[0]
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
        from cStringIO import StringIO
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
        dendropy.treesplit.encode_splits(first, lowest_relevant_bit=2)
        dendropy.treesplit.encode_splits(second, lowest_relevant_bit=2)
        self.assertEquals(2, first.symmetric_difference(second))


if __name__ == "__main__":
    unittest.main()


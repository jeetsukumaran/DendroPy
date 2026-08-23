#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

import unittest
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap
import dendropy

class TreeArrayBasicTreeAccession(unittest.TestCase):

    def get_trees(self, taxon_namespace=None):
        trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path(
                "pythonidae.reference-trees.nexus"),
                "nexus",
                taxon_namespace=taxon_namespace)
        return trees

    def verify_tree_array(self, tree_array, source_trees, ignore_edges=False):
        self.assertEqual(len(tree_array), len(source_trees))
        for idx, source_tree in enumerate(source_trees):
            source_splits = [b.split_bitmask for b in source_tree.encode_bipartitions()]
            tss_splits, tss_edges = tree_array.get_split_bitmask_and_edge_tuple(idx)
            self.assertEqual(len(tss_splits), len(source_splits))
            self.assertEqual(set(tss_splits), set(source_splits))
            # since encoding is done in postorder, we can rely on correspondence of index ...
            for idx, nd in enumerate(source_tree.postorder_node_iter()):
                source_split = nd.edge.bipartition.split_bitmask
                tss_split = tss_splits[idx]
                tss_edge_length = tss_edges[idx]
                self.assertEqual(source_split, tss_split)
                if nd.edge.length is None:
                    self.assertEqual(tss_edge_length, 0)
                else:
                    self.assertAlmostEqual(tss_edge_length, nd.edge.length)

    def test_add_tree(self):
        trees = self.get_trees()
        tree_array = dendropy.TreeArray(taxon_namespace=trees.taxon_namespace)
        for tree in trees:
            tree_array.add_tree(tree)
        self.verify_tree_array(tree_array, trees)


if __name__ == "__main__":
    unittest.main()

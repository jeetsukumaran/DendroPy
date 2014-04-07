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
Tests basic Tree copying etc.
"""

import unittest
import dendropy
import copy
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import compare_and_validate

class TestTreeIdentity(unittest.TestCase):

    def setUp(self):
        self.t1 = dendropy.Tree()
        self.t2 = dendropy.Tree()

    def test_equal(self):
        self.assertNotEqual(self.t1, self.t2)

    def test_hash_dict_membership(self):
        k = {}
        k[self.t1] = 1
        k[self.t2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.t1], 1)
        self.assertEqual(k[self.t2], 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.t1)
        k.add(self.t2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)

class TestTreeCopying(datagen_curated_test_tree.CuratedTestTree, compare_and_validate.Comparator, unittest.TestCase):

    def test_copy(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree()
        for tree2 in (
                tree1.clone(0),
                copy.copy(tree1),
                tree1.clone(1),
                tree1.taxon_namespace_scoped_copy(),
                # dendropy.Tree(tree),
                ):
            self.compare_distinct_trees(tree1, tree2)


if __name__ == "__main__":
    unittest.main()

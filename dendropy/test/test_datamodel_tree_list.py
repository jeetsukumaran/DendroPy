# !/usr/bin/env python

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
Tests for TreeList.
"""

import unittest
from dendropy import TaxonNamespace, Tree, TreeList
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import compare_and_validate

class TreeListIdentity(unittest.TestCase):

    def setUp(self):
        self.tns = TaxonNamespace()
        self.t1 = TreeList(label="a", taxon_namespace=self.tns)
        self.t2 = TreeList(label="a", taxon_namespace=self.tns)
        self.t3 = TreeList(label="a")

    def test_equal(self):
        # two distinct :class:`TreeList` objects are equal
        # if they have the same namespace and trees
        trees = [Tree() for i in range(5)]
        for tree in trees:
            self.t1._trees.append(tree)
            self.t2._trees.append(tree)
        self.assertEqual(self.t1, self.t2)

    def test_unequal1(self):
        # two distinct :class:`TreeList` objects are equal
        # if they have the same namespace and trees
        trees1 = [Tree() for i in range(5)]
        for tree in trees1:
            self.t1._trees.append(tree)
        trees2 = [Tree() for i in range(5)]
        for tree in trees2:
            self.t2._trees.append(tree)
        self.assertNotEqual(self.t1, self.t2)

    def test_unequal2(self):
        # two distinct :class:`TreeList` objects are equal
        # if they have the same namespace and trees
        trees1 = [Tree() for i in range(5)]
        for tree in trees1:
            self.t1._trees.append(tree)
            self.t3._trees.append(tree)
        self.assertNotEqual(self.t1, self.t3)

    def test_hash_dict_membership(self):
        k = {}
        k[self.t1] = 1
        k[self.t2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.t1], 1)
        self.assertEqual(k[self.t2], 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        del k[self.t1]
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.t1)
        k.add(self.t2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        k.discard(self.t1)
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

class TreeListCreatingAndCloning(
        compare_and_validate.Comparator,
        unittest.TestCase):

    def create_with_taxon_namespace(self):
        tns = dendropy.TaxonNamespace()
        tt = TreeList(label="a", taxon_namespace=tns)
        self.assertEqual(tt.label, "a")
        self.assertIs(tt.taxon_namespace, tns)


if __name__ == "__main__":
    unittest.main()

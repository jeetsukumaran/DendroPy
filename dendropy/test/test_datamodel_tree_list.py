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

import sys
import unittest
import dendropy
from dendropy import Taxon, TaxonNamespace, Tree, TreeList
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

class TestTreeListUpdateTaxonNamespace(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        trees = []
        for idx in range(5):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)
            trees.append(tree1)
        self.expected_labels = set()
        self.expected_taxa = set()
        node_label_to_taxon_label_map = {
            "a" : "z01",
            "b" : "<NONE>",
            "c" : "z03",
            "e" : "z04",
            "f" : "z05",
            "g" : "z06",
            "h" : None,
            "i" : None,
            "j" : "z09",
            "k" : "z10",
            "l" : "z11",
            "m" : "<NONE>",
            "n" : None,
            "o" : "z14",
            "p" : "z15",
                }
        registry = {}
        for tree_idx, tree in enumerate(trees):
            for nd in tree:
                if nd.label is not None:
                    if tree_idx > 3:
                        nd.label = node_label_to_taxon_label_map[nd.label]
                    if nd.label == "<NONE>":
                        try:
                            t = registry[None]
                        except KeyError:
                            t = dendropy.Taxon(label=None)
                            registry[None] = t
                        self.expected_labels.add(None)
                    else:
                        try:
                            t = registry[nd.label]
                        except KeyError:
                            t = dendropy.Taxon(label=nd.label)
                            registry[nd.label] = t
                        self.expected_labels.add(nd.label)
                    nd.taxon = t
                    self.expected_taxa.add(nd.taxon)
        self.tree_list = dendropy.TreeList()
        self.tree_list._trees = trees

    def test_noop_update_with_no_taxa(self):
        trees = []
        tns = dendropy.TaxonNamespace()
        for idx in range(5):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True,
                    taxon_namespace=tns)
            trees.append(tree1)
        tlst = dendropy.TreeList(taxon_namespace=tns)
        tlst._trees = trees
        original_tns = tlst.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        tlst.update_taxon_namespace()
        self.assertIs(tlst.taxon_namespace, original_tns)
        for tree in tlst:
            self.assertIs(tree.taxon_namespace, tlst.taxon_namespace)
        self.assertEqual(len(original_tns), 0)

    def test_update(self):
        original_tns = self.tree_list.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        self.tree_list.update_taxon_namespace()
        self.tree_list.update_taxon_namespace()
        self.tree_list.update_taxon_namespace()
        for tree in self.tree_list:
            self.assertIs(tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertIs(self.tree_list.taxon_namespace, original_tns)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        if sys.hexversion < 0x03000000:
            self.assertItemsEqual(new_taxa, self.expected_taxa)
            self.assertItemsEqual(new_labels, self.expected_labels)
        else:
            self.assertCountEqual(new_taxa, self.expected_taxa)
            self.assertCountEqual(new_labels, self.expected_labels)

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

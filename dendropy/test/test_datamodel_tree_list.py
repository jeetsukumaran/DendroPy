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
import collections
import dendropy
import random
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

class TestTreeListMigrateAndReconstructTaxonNamespace(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(8):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True,
                    taxon_namespace=tns)
            trees.append(tree)
        self.node_label_to_taxon_label_map = {
            "a" : "a",
            "b" : "a",
            "c" : "2",
            "e" : "2",
            "f" : "b",
            "g" : "B",
            "h" : "B",
            "i" : "h",
            "j" : "H",
            "k" : "h",
            "l" : None,
            "m" : None,
            "n" : "H",
            "o" : "J",
            "p" : "j",
                }
        self.original_taxa = []
        registry = {}
        for tree in trees:
            for idx, nd in enumerate(tree):
                try:
                    t = registry[nd.label]
                except KeyError:
                    taxon_label = self.node_label_to_taxon_label_map[nd.label]
                    t = dendropy.Taxon(label=taxon_label)
                    registry[nd.label] = t
                    self.original_taxa.append(t)
                tree.taxon_namespace.add_taxon(t)
                nd.taxon = t
                nd.original_taxon = t
        assert len(tree.taxon_namespace) == len(self.node_label_to_taxon_label_map)
        assert len(tree.taxon_namespace) == len(self.original_taxa)
        self.tree_list = dendropy.TreeList(taxon_namespace=tns)
        self.tree_list._trees = trees

    def verify_taxon_namespace_reconstruction(self,
            unify_taxa_by_label=False,
            case_insensitive_label_mapping=False,
            original_tns=None,
            redundant_taxa=False):
        if unify_taxa_by_label:
            if case_insensitive_label_mapping:
                expected_labels = list(set((label.upper() if label is not None else None) for label in self.node_label_to_taxon_label_map.values()))
            else:
                expected_labels = list(set(label for label in self.node_label_to_taxon_label_map.values()))
        else:
            expected_labels = [label for label in self.node_label_to_taxon_label_map.values()]
        for tree in self.tree_list:
            seen_taxa = []
            self.assertIs(tree.taxon_namespace, self.tree_list.taxon_namespace)
            for nd in tree:
                self.assertIsNot(nd.taxon, nd.original_taxon)
                if case_insensitive_label_mapping and nd.taxon.label is not None:
                    self.assertEqual(nd.taxon.label.upper(), nd.original_taxon.label.upper())
                    self.assertEqual(self.node_label_to_taxon_label_map[nd.label].upper(), nd.taxon.label.upper())
                else:
                    self.assertEqual(nd.taxon.label, nd.original_taxon.label)
                    self.assertEqual(self.node_label_to_taxon_label_map[nd.label], nd.taxon.label)
                self.assertNotIn(nd.original_taxon, tree.taxon_namespace)
                self.assertIn(nd.original_taxon, self.original_taxa)
                self.assertIn(nd.taxon, tree.taxon_namespace)
                self.assertNotIn(nd.taxon, self.original_taxa)
                if original_tns is not None:
                    self.assertNotIn(nd.taxon, original_tns)
                if nd.taxon not in seen_taxa:
                    seen_taxa.append(nd.taxon)
                else:
                    self.assertTrue(unify_taxa_by_label or redundant_taxa)
                    if case_insensitive_label_mapping:
                        self.assertIn(nd.taxon.label, [t.label for t in seen_taxa])
                    else:
                        if nd.taxon.label is None:
                            self.assertIs(nd.original_taxon.label, None)
                            self.assertEqual([t.label for t in seen_taxa].count(None), 1)
                        else:
                            x1 = [t.label.upper() for t in seen_taxa if t.label is not None]
                            self.assertIn(nd.taxon.label.upper(), x1)
            self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
            if case_insensitive_label_mapping:
                seen_labels = [(t.label.upper() if t.label is not None else None) for t in seen_taxa]
            else:
                seen_labels = [t.label for t in seen_taxa]
            c1 = collections.Counter(expected_labels)
            c2 = collections.Counter(seen_labels)
            self.assertEqual(c1, c2)

        self.assertEqual(len(tree.taxon_namespace), len(expected_labels))
        if not unify_taxa_by_label and not redundant_taxa:
            self.assertEqual(len(tree.taxon_namespace), len(self.node_label_to_taxon_label_map))

    def test_basic_reconstruction(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(5):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False,
                    taxon_namespace=tns)
            trees.append(tree)
        tree_list = dendropy.TreeList(taxon_namespace=tns)
        tree_list._trees = trees
        new_tns = dendropy.TaxonNamespace()
        tree_list.taxon_namespace = new_tns
        tree_list.reconstruct_taxon_namespace(
                unify_taxa_by_label=False,
                case_insensitive_label_mapping=False)
        self.assertIsNot(tree_list.taxon_namespace, tns)
        self.assertIs(tree_list.taxon_namespace, new_tns)
        self.assertEqual(len(tree_list.taxon_namespace), len(tns))
        original_labels = [t.label for t in tns]
        new_labels = [t.label for t in new_tns]
        if sys.hexversion < 0x03000000:
            self.assertItemsEqual(new_labels, original_labels)
        else:
            self.assertCountEqual(new_labels, original_labels)
        for tree in tree_list:
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            for nd in tree:
                if nd.taxon is not None:
                    self.assertIn(nd.taxon, tree.taxon_namespace)
                    self.assertNotIn(nd.taxon, tns)

    def test_reconstruct_taxon_namespace_non_unifying(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=False,
                case_insensitive_label_mapping=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_insensitive_label_mapping=False)

    def test_reconstruct_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True,
                case_insensitive_label_mapping=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=False,
                original_tns=original_tns)

    def test_reconstruct_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True,
                case_insensitive_label_mapping=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=True,
                original_tns=original_tns)

    def test_basic_migration(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(5):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False,
                    taxon_namespace=tns)
            trees.append(tree)
        tree_list = dendropy.TreeList(taxon_namespace=tns)
        tree_list._trees = trees
        new_tns = dendropy.TaxonNamespace()
        tree_list.taxon_namespace = new_tns
        tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False,
                case_insensitive_label_mapping=False)
        self.assertIsNot(tree_list.taxon_namespace, tns)
        self.assertIs(tree_list.taxon_namespace, new_tns)
        self.assertEqual(len(tree_list.taxon_namespace), len(tns))
        original_labels = [t.label for t in tns]
        new_labels = [t.label for t in new_tns]
        if sys.hexversion < 0x03000000:
            self.assertItemsEqual(new_labels, original_labels)
        else:
            self.assertCountEqual(new_labels, original_labels)
        for tree in tree_list:
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            for nd in tree:
                if nd.taxon is not None:
                    self.assertIn(nd.taxon, tree.taxon_namespace)
                    self.assertNotIn(nd.taxon, tns)

    def test_migrate_taxon_namespace_non_unifying(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False,
                case_insensitive_label_mapping=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_insensitive_label_mapping=False,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=False,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_insensitive_label_mapping=True,
                original_tns=original_tns)

class TestTreeListAppend(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        self.native_tns = dendropy.TaxonNamespace()
        self.tree_list = dendropy.TreeList(taxon_namespace=self.native_tns)
        self.foreign_tns = dendropy.TaxonNamespace()
        self.foreign_tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                taxon_namespace=self.foreign_tns)
        for nd in self.foreign_tree:
            nd.original_taxon = nd.taxon
        self.check_tns = dendropy.TaxonNamespace()
        self.check_tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                taxon_namespace=self.check_tns)

    def test_append_default(self):
        self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
        self.tree_list.append(self.foreign_tree)
        self.assertEqual(len(self.tree_list), 1)
        self.assertIn(self.foreign_tree, self.tree_list)
        self.assertIs(self.foreign_tree, self.tree_list[0])
        self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
        self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
        for nd in self.foreign_tree:
            if nd.taxon:
                self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
                self.assertIsNot(nd.taxon, nd.original_taxon)
                self.assertIn(nd.original_taxon, self.foreign_tns)
                self.assertNotIn(nd.original_taxon, self.tree_list.taxon_namespace)
                self.assertEqual(nd.taxon.label, nd.original_taxon.label)

    def test_append_add(self):
        self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
        self.tree_list.append(self.foreign_tree, taxon_import_strategy="update")
        self.assertEqual(len(self.tree_list), 1)
        self.assertIn(self.foreign_tree, self.tree_list)
        self.assertIs(self.foreign_tree, self.tree_list[0])
        self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
        self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
        for nd in self.foreign_tree:
            if nd.taxon:
                self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
                self.assertIs(nd.taxon, nd.original_taxon)
                self.assertIn(nd.original_taxon, self.foreign_tns)
                self.assertIn(nd.original_taxon, self.tree_list.taxon_namespace)

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

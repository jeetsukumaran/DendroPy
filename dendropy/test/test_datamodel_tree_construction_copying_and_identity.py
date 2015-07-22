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
Tests basic Tree copying etc.
"""

import unittest
import dendropy
import copy
from dendropy.test.support import curated_test_tree
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

class TestTreeCopying(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.Comparator,
        unittest.TestCase):

    def add_annotations(self, tree):
        for idx, nd in enumerate(tree):
            if idx % 2 == 0:
                nd.edge.label = "E{}".format(idx)
                nd.edge.length = idx
            an1 = nd.annotations.add_new("a{}".format(idx),
                    "{}{}{}".format(nd.label, nd.taxon, idx))
            an2 = nd.annotations.add_bound_attribute("label")
            an3 = an1.annotations.add_bound_attribute("name")
            ae1 = nd.edge.annotations.add_new("a{}".format(idx),
                    "{}{}".format(nd.edge.label, idx))
            ae2 = nd.edge.annotations.add_bound_attribute("label")
            ae3 = ae1.annotations.add_bound_attribute("name")
        tree.annotations.add_new("a", 0)
        tree.label = "hello"
        b = tree.annotations.add_bound_attribute("label")
        b.annotations.add_new("c", 3)
        for idx, taxon in enumerate(tree.taxon_namespace):
            a = taxon.annotations.add_new("!color", str(idx))
            a.annotations.add_new("setbytest", "a")

    def test_copy(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.add_annotations(tree1)
        for tree2 in (
                # tree1.clone(0),
                # copy.copy(tree1),
                # tree1.clone(1),
                # tree1.taxon_namespace_scoped_copy(),
                dendropy.Tree(tree1),
                ):
            self.compare_distinct_trees(tree1, tree2,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=True,
                    compare_taxon_annotations=False)
            # Redundant, given the above
            # But for sanity's sake ...
            nodes1 = [nd for nd in tree1]
            nodes2 = [nd for nd in tree2]
            self.assertEqual(len(nodes1), len(nodes2))
            for nd1, nd2 in zip(nodes1, nodes2):
                self.assertIsNot(nd1, nd2)
                self.assertEqual(nd1.label, nd2.label)
                self.assertIs(nd1.taxon, nd2.taxon)

    def test_deepcopy_including_namespace(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.add_annotations(tree1)
        for idx, tree2 in enumerate((
                tree1.clone(2),
                copy.deepcopy(tree1),
                )):
            self.compare_distinct_trees(tree1, tree2,
                    taxon_namespace_scoped=False,
                    compare_tree_annotations=True,
                    compare_taxon_annotations=False)
            # Redundant, given the above
            # But for sanity's sake ...
            nodes1 = [nd for nd in tree1]
            nodes2 = [nd for nd in tree2]
            self.assertEqual(len(nodes1), len(nodes2))
            for nd1, nd2 in zip(nodes1, nodes2):
                self.assertIsNot(nd1, nd2)
                self.assertEqual(nd1.label, nd2.label)
                self.assertIsNot(nd1.taxon, nd2.taxon)
                self.assertEqual(nd1.taxon.label, nd2.taxon.label)

    def test_deepcopy_excluding_namespace(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.add_annotations(tree1)
        tree2 = dendropy.Tree(tree1, taxon_namespace=dendropy.TaxonNamespace())
        self.compare_distinct_trees(tree1, tree2,
                taxon_namespace_scoped=False,
                compare_tree_annotations=True,
                compare_taxon_annotations=False)

class TestSpecialTreeConstruction(
        curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def test_construction_from_another_tree_different_label(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        tree1.label = "tree1"
        self.assertEqual(tree1.label, "tree1")
        tree2 = dendropy.Tree(tree1, label="tree2")
        self.assertEqual(tree2.label, "tree2")
        self.assertNotEqual(tree1.label, "tree2")
        self.assertNotEqual(tree1.label, tree2.label)

    def test_construction_from_given_seed_node(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        tree2 = dendropy.Tree(seed_node=tree1.seed_node)
        self.assertIs(tree2.seed_node, tree1.seed_node)

    def test_construction_from_given_seed_node(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        tree2 = dendropy.Tree(seed_node=tree1.seed_node)
        self.assertIs(tree2.seed_node, tree1.seed_node)
        self.assertIsNot(tree1.taxon_namespace, tree2.taxon_namespace)
        self.assertEqual(len(tree1.taxon_namespace), len(tree2.taxon_namespace))
        for taxon1 in tree1.taxon_namespace:
            self.assertIn(taxon1, tree2.taxon_namespace)
        for taxon2 in tree2.taxon_namespace:
            self.assertIn(taxon2, tree1.taxon_namespace)
        for nd in tree2:
            self.assertIn(nd.taxon, tree2.taxon_namespace)

    def test_cloning_construction_with_taxon_namespace(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        tns = dendropy.TaxonNamespace()
        tree2 = dendropy.Tree(tree1, taxon_namespace=tns)
        self.assertIs(tree2.taxon_namespace, tns)
        self.assertIsNot(tree1.taxon_namespace, tree2.taxon_namespace)
        self.assertEqual(len(tree1.taxon_namespace), len(tree2.taxon_namespace))
        for nd in tree2:
            self.assertIn(nd.taxon, tree2.taxon_namespace)
            self.assertNotIn(nd.taxon, tree1.taxon_namespace)
        for nd in tree1:
            self.assertIn(nd.taxon, tree1.taxon_namespace)
            self.assertNotIn(nd.taxon, tree2.taxon_namespace)

if __name__ == "__main__":
    unittest.main()

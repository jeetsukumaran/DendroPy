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
Tests Tree taxon management
"""

import unittest
import dendropy
import collections
import copy
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import compare_and_validate

class TestTreeUpdateTaxonNamespace(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        self.tree1, self.anodes1, self.lnodes1, self.inodes1 = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True)
        self.expected_labels = set([nd.label for nd in self.anodes1 if nd.label is not None])
        self.expected_taxa = set()
        for nd in self.tree1:
            if nd.label is not None:
                nd.taxon = dendropy.Taxon(label=nd.label)
                self.expected_taxa.add(nd.taxon)
        assert len(self.expected_labels) == len(self.anodes1)
        assert len(self.expected_taxa) == len(self.expected_labels)

    # def test_infer_taxa(self):
    #     tns1 = self.tree1.taxon_namespace
    #     self.assertEqual(len(tns1), 0)
    #     tns2 = self.tree1.infer_taxa()
    #     self.assertNotEqual(len(tns2), 0)
    #     self.assertIs(self.tree1.taxon_namespace, tns2)
    #     self.assertIsNot(self.tree1.taxon_namespace, tns1)
    #     self.assertEqual(set(tns2._taxa), self.expected_taxa)

    # def test_reindex_subcomponent_taxa(self):
    #     tns1 = self.tree1.taxon_namespace
    #     self.assertEqual(len(tns1), 0)
    #     self.tree1.reindex_subcomponent_taxa()
    #     tns2 = self.tree1.taxon_namespace
    #     # self.assertIsNot(tns1, tns2) tns1 *IS* tns2
    #     self.assertEqual(len(tns2), len(self.expected_taxa))
    #     self.assertEqual(set(t.label for t in self.tree1.taxon_namespace),
    #             self.expected_labels)

    def test_update_taxon_namespace(self):
        tns1 = self.tree1.taxon_namespace
        self.assertEqual(len(tns1), 0)
        tns2 = self.tree1.update_taxon_namespace()
        self.assertIs(self.tree1.taxon_namespace, tns1)
        self.assertEqual(set(tns2._taxa), self.expected_taxa)
        self.assertEqual(len(tns2._taxa), len(self.expected_labels))

class TestTreeReconstructTaxonNamespace(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        self.tree, self.anodes, self.lnodes, self.inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True)
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
        self.label_to_taxa = collections.defaultdict(list)
        for nd in self.tree:
            taxon_label = self.node_label_to_taxon_label_map[nd.label]
            t = dendropy.Taxon(label=taxon_label)
            self.tree.taxon_namespace.add_taxon(t)
            nd.taxon = t
            nd.original_taxon = t
            self.original_taxa.append(t)
            self.label_to_taxa[t.label].append(t)
        assert len(self.tree.taxon_namespace) == len(self.node_label_to_taxon_label_map)
        assert len(self.tree.taxon_namespace) == len(self.original_taxa)

    def test_reconstruct_taxon_namespace_non_unifying(self):
        original_tns = self.tree.taxon_namespace
        self.tree.taxon_namespace = dendropy.TaxonNamespace()
        self.assertEqual(len(self.tree.taxon_namespace), 0)
        self.tree.reconstruct_taxon_namespace(unify_taxa_by_label=False)
        seen_taxa = []
        for nd in self.tree:
            seen_taxa.append(nd.taxon)
            self.assertIsNot(nd.taxon, nd.original_taxon)
            self.assertEqual(nd.taxon.label, nd.original_taxon.label)
            self.assertNotIn(nd.original_taxon, self.tree.taxon_namespace)
            self.assertIn(nd.original_taxon, self.original_taxa)
            self.assertIn(nd.taxon, self.tree.taxon_namespace)
            self.assertNotIn(nd.taxon, self.original_taxa)
            self.assertEqual(self.node_label_to_taxon_label_map[nd.label], nd.taxon.label)
        self.assertEqual(len(self.tree.taxon_namespace), len(self.node_label_to_taxon_label_map))

    def test_reconstruct_taxon_namespace_unifying_case_sensitive(self):
        pass

    def test_reconstruct_taxon_namespace_unifying_case_insensitive(self):
        pass

    def test_reconstruct_taxon_namespace_with_taxon_mapping(self):
        pass

if __name__ == "__main__":
    unittest.main()

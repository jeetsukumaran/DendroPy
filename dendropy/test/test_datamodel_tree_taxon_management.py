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
Tests Tree taxon management
"""

import os
import unittest
import dendropy
import collections
import copy
from dendropy.test.support import curated_test_tree
from dendropy.test.support import compare_and_validate
from dendropy.test.support import dendropytest

class TestTreeUpdateTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree1, self.anodes1, self.lnodes1, self.inodes1 = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
        self.expected_labels = set([nd.label for nd in self.anodes1 if nd.label is not None])
        self.expected_taxa = set()
        for nd in self.tree1:
            if nd.label is not None:
                nd.taxon = dendropy.Taxon(label=nd.label)
                self.expected_taxa.add(nd.taxon)
        assert len(self.expected_labels) == len(self.anodes1)
        assert len(self.expected_taxa) == len(self.expected_labels)

    def test_noop_update_with_no_taxa(self):
        tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
        original_tns = tree.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        tree.update_taxon_namespace()
        self.assertIs(tree.taxon_namespace, original_tns)
        self.assertEqual(len(original_tns), 0)

    def test_noop_update(self):
        tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        original_tns = tree.taxon_namespace
        original_taxa = [t for t in original_tns]
        original_labels = [t.label for t in original_tns]
        tree.update_taxon_namespace()
        self.assertIs(tree.taxon_namespace, original_tns)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        self.assertEqual(new_taxa, original_taxa)
        self.assertEqual(new_labels, original_labels)

    def test_update_taxon_namespace(self):
        tns1 = self.tree1.taxon_namespace
        self.assertEqual(len(tns1), 0)
        tns2 = self.tree1.update_taxon_namespace()
        self.assertIs(self.tree1.taxon_namespace, tns1)
        self.assertEqual(set(tns2._taxa), self.expected_taxa)
        self.assertEqual(len(tns2._taxa), len(self.expected_labels))

class TestTreeMigrateAndReconstructTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree, self.anodes, self.lnodes, self.inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
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
        for idx, nd in enumerate(self.tree):
            taxon_label = self.node_label_to_taxon_label_map[nd.label]
            t = dendropy.Taxon(label=taxon_label)
            self.tree.taxon_namespace.add_taxon(t)
            nd.taxon = t
            nd.original_taxon = t
            self.original_taxa.append(t)
        assert len(self.tree.taxon_namespace) == len(self.node_label_to_taxon_label_map)
        assert len(self.tree.taxon_namespace) == len(self.original_taxa)

    def create_redundant_taxa(self):
        nodes_to_unify_taxa = []
        for nd in self.tree:
            if nd.taxon.label and nd.taxon.label.upper() == "J":
                nodes_to_unify_taxa.append(nd)
        assert len(nodes_to_unify_taxa) >= 2
        utaxon = nodes_to_unify_taxa[0].taxon
        utaxon.node_label = nodes_to_unify_taxa[0].label
        for nd in nodes_to_unify_taxa[1:]:
            self.tree.taxon_namespace.remove_taxon(nd.taxon)
            self.original_taxa.remove(nd.taxon)
            nd.original_taxon = utaxon
            nd.taxon = utaxon
            del self.node_label_to_taxon_label_map[nd.label]
            nd.label = utaxon.node_label

    def verify_taxon_namespace_reconstruction(self,
            unify_taxa_by_label=False,
            case_sensitive_label_mapping=True,
            original_tns=None,
            redundant_taxa=False):
        seen_taxa = []
        if unify_taxa_by_label:
            if not case_sensitive_label_mapping:
                expected_labels = []
                for label in self.node_label_to_taxon_label_map.values():
                    if label is None:
                        expected_labels.append(label)
                    else:
                        label = label.upper()
                        if label not in expected_labels:
                            expected_labels.append(label)
            else:
                expected_labels = list(set(label for label in self.node_label_to_taxon_label_map.values()))
        else:
            expected_labels = [label for label in self.node_label_to_taxon_label_map.values()]
        for nd in self.tree:
            self.assertIsNot(nd.taxon, nd.original_taxon)
            if (not case_sensitive_label_mapping) and nd.taxon.label is not None:
                self.assertEqual(nd.taxon.label.upper(), nd.original_taxon.label.upper())
                self.assertEqual(self.node_label_to_taxon_label_map[nd.label].upper(), nd.taxon.label.upper())
            else:
                self.assertEqual(nd.taxon.label, nd.original_taxon.label)
                self.assertEqual(self.node_label_to_taxon_label_map[nd.label], nd.taxon.label)
            self.assertNotIn(nd.original_taxon, self.tree.taxon_namespace)
            self.assertIn(nd.original_taxon, self.original_taxa)
            self.assertIn(nd.taxon, self.tree.taxon_namespace)
            self.assertNotIn(nd.taxon, self.original_taxa)
            if original_tns is not None:
                self.assertNotIn(nd.taxon, original_tns)
            if nd.taxon not in seen_taxa:
                seen_taxa.append(nd.taxon)
            else:
                self.assertTrue(unify_taxa_by_label or redundant_taxa)
                if not case_sensitive_label_mapping:
                    self.assertIn(nd.taxon.label, [t.label for t in seen_taxa])
                else:
                    if nd.taxon.label is None:
                        self.assertIs(nd.original_taxon.label, None)
                        self.assertEqual([t.label for t in seen_taxa].count(None), 1)
                    else:
                        x1 = [t.label.upper() for t in seen_taxa if t.label is not None]
                        self.assertIn(nd.taxon.label.upper(), x1)
        self.assertEqual(len(self.tree.taxon_namespace), len(expected_labels))
        if not unify_taxa_by_label and not redundant_taxa:
            self.assertEqual(len(self.tree.taxon_namespace), len(self.node_label_to_taxon_label_map))
        self.assertEqual(len(seen_taxa), len(self.tree.taxon_namespace))
        if not case_sensitive_label_mapping:
            seen_labels = [(t.label.upper() if t.label is not None else None) for t in seen_taxa]
        else:
            seen_labels = [t.label for t in seen_taxa]
        c1 = collections.Counter(expected_labels)
        c2 = collections.Counter(seen_labels)
        self.assertEqual(c1, c2)

    def test_noop_taxon_namespace_reconstruction(self):
        tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        original_tns = tree.taxon_namespace
        original_tns.is_case_sensitive = True
        original_taxa = [t for t in original_tns]
        original_labels = [t.label for t in original_tns]
        tree.reconstruct_taxon_namespace(
                unify_taxa_by_label=False)
        self.assertIs(tree.taxon_namespace, original_tns)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        self.assertEqual(new_taxa, original_taxa)
        self.assertEqual(new_labels, original_labels)

    def test_reconstruct_taxon_namespace_non_unifying(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree._taxon_namespace = new_tns
        self.assertEqual(len(self.tree.taxon_namespace), 0)
        self.tree.reconstruct_taxon_namespace(unify_taxa_by_label=False)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True)

    def test_reconstruct_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree._taxon_namespace = new_tns
        self.assertEqual(len(self.tree.taxon_namespace), 0)
        self.tree.reconstruct_taxon_namespace(unify_taxa_by_label=True)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_reconstruct_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = False
        self.tree._taxon_namespace = new_tns
        self.assertEqual(len(self.tree.taxon_namespace), 0)
        self.tree.reconstruct_taxon_namespace(unify_taxa_by_label=True)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=False,
                original_tns=original_tns)

    def test_reconstruct_taxon_namespace_with_redundant_taxa(self):
        for (unify, ci) in [
                (False, True),
                (True, True),
                (True, False), ]:
            self.setUp()
            self.create_redundant_taxa()
            original_tns = self.tree.taxon_namespace
            new_tns = dendropy.TaxonNamespace()
            new_tns.is_case_sensitive = ci
            self.tree._taxon_namespace = new_tns
            self.assertEqual(len(self.tree.taxon_namespace), 0)
            self.tree.reconstruct_taxon_namespace(unify_taxa_by_label=unify)
            self.assertIsNot(self.tree.taxon_namespace, original_tns)
            self.assertIs(self.tree.taxon_namespace, new_tns)
            self.verify_taxon_namespace_reconstruction(
                    unify_taxa_by_label=unify,
                    case_sensitive_label_mapping=ci,
                    original_tns=original_tns,
                    redundant_taxa=True)

    def test_reconstruct_taxon_namespace_mapping(self):
        for (unify, ci) in [
                (False, False),
                (True, False),
                (True, True), ]:
            self.setUp()
            original_tns = self.tree.taxon_namespace
            new_tns = dendropy.TaxonNamespace()
            new_tns.is_case_sensitive = ci
            self.tree._taxon_namespace = new_tns
            self.assertEqual(len(self.tree.taxon_namespace), 0)
            memo = {}
            for taxon in self.original_taxa:
                memo[taxon] = dendropy.Taxon()
            memo_copy = dict(memo)
            self.tree.reconstruct_taxon_namespace(
                    unify_taxa_by_label=unify,
                    taxon_mapping_memo=memo)
            self.assertIsNot(self.tree.taxon_namespace, original_tns)
            self.assertIs(self.tree.taxon_namespace, new_tns)
            for nd in self.tree:
                self.assertIs(nd.taxon, memo_copy[nd.original_taxon])

    def test_noop_migrate_taxon_namespace(self):
        tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        original_tns = tree.taxon_namespace
        original_tns.is_case_sensitive = True
        original_taxa = [t for t in original_tns]
        original_labels = [t.label for t in original_tns]
        tree.migrate_taxon_namespace(
                original_tns,
                unify_taxa_by_label=False)
        self.assertIs(tree.taxon_namespace, original_tns)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        self.assertEqual(new_taxa, original_taxa)
        self.assertEqual(new_labels, original_labels)

    def test_simple_migrate_taxon_namespace(self):
        tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
        original_tns = tree.taxon_namespace
        original_taxa = [t for t in original_tns]
        original_labels = [t.label for t in original_tns]
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        tree.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False)
        self.assertIsNot(tree.taxon_namespace, original_tns)
        self.assertIs(tree.taxon_namespace, new_tns)
        new_taxa = [t for t in new_tns]
        self.assertEqual(len(new_taxa), len(original_taxa))
        for t1, t2 in zip(new_taxa, original_taxa):
            self.assertIsNot(t1, t2)
            self.assertEqual(t1.label, t2.label)

    def test_migrate_taxon_namespace_non_unifying(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = False
        self.tree.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True)
        self.assertIsNot(self.tree.taxon_namespace, original_tns)
        self.assertIs(self.tree.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=False,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_mapping(self):
        for (unify, ci) in [
                (False, True),
                (True, True),
                (True, False), ]:
            self.setUp()
            original_tns = self.tree.taxon_namespace
            new_tns = dendropy.TaxonNamespace()
            new_tns.is_case_sensitive = ci
            memo = {}
            for taxon in self.original_taxa:
                memo[taxon] = dendropy.Taxon()
            memo_copy = dict(memo)
            self.tree.migrate_taxon_namespace(
                    new_tns,
                    unify_taxa_by_label=unify,
                    taxon_mapping_memo=memo)
            self.assertIsNot(self.tree.taxon_namespace, original_tns)
            self.assertIs(self.tree.taxon_namespace, new_tns)
            for nd in self.tree:
                self.assertIs(nd.taxon, memo_copy[nd.original_taxon])

    def test_migrate_taxon_namespace_with_redundant_taxa(self):
        for (unify, ci) in [
                (False, True),
                (True, True),
                (True, False), ]:
            self.setUp()
            self.create_redundant_taxa()
            original_tns = self.tree.taxon_namespace
            new_tns = dendropy.TaxonNamespace()
            new_tns.is_case_sensitive = ci
            self.tree.migrate_taxon_namespace(
                    new_tns,
                    unify_taxa_by_label=unify)
            self.assertIsNot(self.tree.taxon_namespace, original_tns)
            self.assertIs(self.tree.taxon_namespace, new_tns)
            self.verify_taxon_namespace_reconstruction(
                    unify_taxa_by_label=unify,
                    case_sensitive_label_mapping=ci,
                    original_tns=original_tns,
                    redundant_taxa=True)

    def test_unassign_taxa(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_randomly_assign_taxa(self):
        self.assertFalse(self.fail_incomplete_tests())

class TestTreeTaxa(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree1, self.anodes1, self.lnodes1, self.inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.expected_taxa = set([nd.taxon for nd in self.anodes1 if nd.taxon is not None])

    def test_basic_taxa(self):
        self.assertEqual(self.tree1.poll_taxa(), self.expected_taxa)

class TestTreePurgeTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree1, self.anodes1, self.lnodes1, self.inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.expected_taxa = set([nd.taxon for nd in self.anodes1 if nd.taxon is not None])

    def test_noop_purge(self):
        self.assertEqual(set(self.tree1.taxon_namespace), self.expected_taxa)
        self.tree1.purge_taxon_namespace()
        self.assertEqual(set(self.tree1.taxon_namespace), self.expected_taxa)

    def test_basic_purge(self):
        self.assertEqual(set(self.tree1.taxon_namespace), self.expected_taxa)
        added_taxa = set(self.expected_taxa)
        for label in ("z1", "z2", "z3", "z4"):
            t = self.tree1.taxon_namespace.new_taxon(label=label)
            added_taxa.add(t)
        self.assertEqual(set(self.tree1.taxon_namespace), added_taxa)
        self.tree1.purge_taxon_namespace()
        self.assertEqual(set(self.tree1.taxon_namespace), self.expected_taxa)

if __name__ == "__main__":
    unittest.main()

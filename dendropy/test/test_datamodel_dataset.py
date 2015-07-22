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
Tests basic DataSet curation
"""

import collections
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import curated_test_tree_list
from dendropy.test.support import standard_file_test_chars

class EmptyDataSetWithLabelCreationTestCase(dendropytest.ExtendedTestCase):

    def test_basic_create_with_label(self):
        ds = dendropy.DataSet(label="d1")
        self.assertEqual(ds.label, "d1")

class DataSetAddTestCase(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_chars.DnaTestChecker.build()
        standard_file_test_chars.ProteinTestChecker.build()

    def setUp(self):
        self.expected_taxon_namespaces = []
        self.standalone_taxon_namespaces = []
        self.standalone_taxon_namespaces.append(dendropy.TaxonNamespace(["t1", "t2", "t3"]))
        self.standalone_taxon_namespaces.append(dendropy.TaxonNamespace(["s1", "s2", "s3"]))
        self.expected_taxon_namespaces.extend(self.standalone_taxon_namespaces)
        self.expected_tree_lists = collections.OrderedDict()
        for i in range(2):
            pdo1 = curated_test_tree_list.get_tree_list(4)
            self.expected_tree_lists[pdo1] = pdo1.taxon_namespace
            self.expected_taxon_namespaces.append(pdo1.taxon_namespace)
            for j in range(2):
                pdo2 = curated_test_tree_list.get_tree_list(4,
                        taxon_namespace=pdo1.taxon_namespace)
                self.expected_tree_lists[pdo2] = pdo2.taxon_namespace

        self.expected_char_matrices = collections.OrderedDict()
        for i in range(2):
            pdo1 = standard_file_test_chars.DnaTestChecker.get_char_matrix_from_class_data()
            self.expected_char_matrices[pdo1] = pdo1.taxon_namespace
            self.expected_taxon_namespaces.append(pdo1.taxon_namespace)
            for j in range(2):
                pdo2 = standard_file_test_chars.ProteinTestChecker.get_char_matrix_from_class_data(taxon_namespace=pdo1.taxon_namespace)
                self.expected_char_matrices[pdo2] = pdo2.taxon_namespace

    def test_basic_add_taxon_namespace(self):
        expected_fundamental_states = set()
        ds = dendropy.DataSet()
        for tns in self.expected_taxon_namespaces:
            ds.add_taxon_namespace(tns)
        self.assertEqual(len(ds.taxon_namespaces), len(self.expected_taxon_namespaces))
        for x1, x2 in zip(ds.taxon_namespaces, self.expected_taxon_namespaces):
            self.assertIs(x1, x2)

    def test_basic_add_tree_list(self):
        ds = dendropy.DataSet()
        expected_taxon_namespaces = collections.OrderedDict()
        for tree_list in self.expected_tree_lists:
            ds.add_tree_list(tree_list)
            expected_taxon_namespaces[self.expected_tree_lists[tree_list]] = True
        self.assertEqual(len(ds.taxon_namespaces), len(expected_taxon_namespaces))
        for x1, x2 in zip(ds.taxon_namespaces, expected_taxon_namespaces):
            self.assertIs(x1, x2)
        self.assertEqual(len(ds.tree_lists), len(self.expected_tree_lists))
        for x1, x2 in zip(ds.tree_lists, self.expected_tree_lists):
            self.assertIs(x1, x2)
            self.assertIs(x1.taxon_namespace, self.expected_tree_lists[x1])
            for t in x1:
                self.assertIs(t.taxon_namespace, x1.taxon_namespace)

    def test_basic_add_char_matrix(self):
        ds = dendropy.DataSet()
        expected_taxon_namespaces = collections.OrderedDict()
        for char_matrix in self.expected_char_matrices:
            ds.add_char_matrix(char_matrix)
            expected_taxon_namespaces[self.expected_char_matrices[char_matrix]] = True
        self.assertEqual(len(ds.taxon_namespaces), len(expected_taxon_namespaces))
        for x1, x2 in zip(ds.taxon_namespaces, expected_taxon_namespaces):
            self.assertIs(x1, x2)
        self.assertEqual(len(ds.char_matrices), len(self.expected_char_matrices))
        for x1, x2 in zip(ds.char_matrices, self.expected_char_matrices):
            self.assertIs(x1, x2)
            self.assertIs(x1.taxon_namespace, self.expected_char_matrices[x1])

    def test_basic_add(self):
        ds = dendropy.DataSet()
        for tns in self.standalone_taxon_namespaces:
            ds.add(tns)
        for tree_list in self.expected_tree_lists:
            ds.add(tree_list)
        for char_matrix in self.expected_char_matrices:
            ds.add(char_matrix)
        self.assertEqual(len(ds.taxon_namespaces), len(self.expected_taxon_namespaces))
        for x1, x2 in zip(ds.taxon_namespaces, self.expected_taxon_namespaces):
            self.assertIs(x1, x2)
        for x1, x2 in zip(ds.tree_lists, self.expected_tree_lists):
            self.assertIs(x1, x2)
        for x1, x2 in zip(ds.char_matrices, self.expected_char_matrices):
            self.assertIs(x1, x2)

    def test_construction(self):
        item_list = []
        item_list.extend(self.standalone_taxon_namespaces)
        item_list.extend(self.expected_tree_lists)
        item_list.extend(self.expected_char_matrices)
        ds = dendropy.DataSet(item_list)
        self.assertEqual(len(ds.taxon_namespaces), len(self.expected_taxon_namespaces))
        for x1, x2 in zip(ds.taxon_namespaces, self.expected_taxon_namespaces):
            self.assertIs(x1, x2)
        for x1, x2 in zip(ds.tree_lists, self.expected_tree_lists):
            self.assertIs(x1, x2)
        for x1, x2 in zip(ds.char_matrices, self.expected_char_matrices):
            self.assertIs(x1, x2)

class DataSetNewTestCase(dendropytest.ExtendedTestCase):

    def test_basic_new_taxon_namespace(self):
        ds = dendropy.DataSet()
        tax_labels = ["a", "b", "c", "d", "e"]
        tns_labels = ["t1", "t2", "t3"]
        tns_list = []
        for tns_label in tns_labels:
             tns = ds.new_taxon_namespace(tax_labels, label=tns_label)
             self.assertTrue(isinstance(tns, dendropy.TaxonNamespace))
             tns_list.append(tns)
        self.assertEqual(len(tns_list), len(tns_labels))
        for tns, tns_label in zip(tns_list, tns_labels):
            self.assertEqual(tns.label, tns_label)
            self.assertEqual(len(tns), len(tax_labels))
            for taxon, taxon_label in zip(tns, tax_labels):
                self.assertEqual(taxon.label, taxon_label)

    def test_basic_new_tree_list(self):
        ds = dendropy.DataSet()
        item_labels = ["a", "b", "c", "d", "e"]
        item_list = []
        for item_idx, item_label in enumerate(item_labels):
            item = ds.new_tree_list(label=item_label)
            item_list.append(item)
        self.assertEqual(len(ds.tree_lists), len(item_labels))
        self.assertEqual(len(ds.tree_lists), len(item_list))
        for t1, t2, label in zip(ds.tree_lists, item_list, item_labels):
            self.assertTrue(isinstance(t1, dendropy.TreeList))
            self.assertIs(t1, t2)
            self.assertEqual(t1.label, label)

    def test_basic_new_char_matrix(self):
        ds = dendropy.DataSet()
        item_labels = ["a", "b", "c", "d", "e", "f"]
        cm_type = [
                "dna",
                "protein",
                "standard",
                dendropy.DnaCharacterMatrix,
                dendropy.ProteinCharacterMatrix,
                dendropy.StandardCharacterMatrix,
                ]
        expected_cm_types = [
                dendropy.DnaCharacterMatrix,
                dendropy.ProteinCharacterMatrix,
                dendropy.StandardCharacterMatrix,
                dendropy.DnaCharacterMatrix,
                dendropy.ProteinCharacterMatrix,
                dendropy.StandardCharacterMatrix,
                ]
        item_list = []
        for item_label, cm_type in zip(item_labels, cm_type):
            item = ds.new_char_matrix(label=item_label,
                    char_matrix_type=cm_type)
            item_list.append(item)
        self.assertEqual(len(ds.char_matrices), len(item_labels))
        self.assertEqual(len(ds.char_matrices), len(item_list))
        for t1, t2, label, expected_cm_types in zip(ds.char_matrices, item_list, item_labels, expected_cm_types):
            self.assertTrue(isinstance(t1, expected_cm_types))
            self.assertIs(t1, t2)
            self.assertEqual(t1.label, label)

class DataSetAttachedTaxonNamespaceModeTestCase(dendropytest.ExtendedTestCase):

    def test_attached_taxon_namespace_default(self):
        ds = dendropy.DataSet()
        tns = dendropy.TaxonNamespace()
        ds.attach_taxon_namespace(tns)
        self.assertTrue(isinstance(tns, dendropy.TaxonNamespace))
        self.assertEqual(len(ds.taxon_namespaces), 1)
        self.assertIn(tns, ds.taxon_namespaces)
        self.assertIs(ds.taxon_namespaces[0], tns)
        self.assertIs(ds.attached_taxon_namespace, tns)
        tns2 = ds.detach_taxon_namespace()
        self.assertIs(ds.attached_taxon_namespace, None)
        self.assertIs(tns2, tns)

    def test_attached_taxon_namespace_new_tree_list(self):
        ds = dendropy.DataSet()
        tns = dendropy.TaxonNamespace()
        ds.attach_taxon_namespace(tns)
        tree_list = ds.new_tree_list(label="q")
        self.assertEqual(tree_list.label, "q")
        self.assertIn(tree_list, ds.tree_lists)
        self.assertIs(tree_list.taxon_namespace, ds.attached_taxon_namespace)
        self.assertEqual(len(ds.taxon_namespaces), 1)

    def test_attached_taxon_namespace_new_char_matrix(self):
        ds = dendropy.DataSet()
        tns = dendropy.TaxonNamespace()
        ds.attach_taxon_namespace(tns)
        char_matrix = ds.new_char_matrix(label="q", char_matrix_type="dna")
        self.assertEqual(char_matrix.label, "q")
        self.assertIn(char_matrix, ds.char_matrices)
        self.assertIs(char_matrix.taxon_namespace, ds.attached_taxon_namespace)
        self.assertEqual(len(ds.taxon_namespaces), 1)

if __name__ == "__main__":
    unittest.main()


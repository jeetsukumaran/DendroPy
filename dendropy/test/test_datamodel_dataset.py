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

class DataSetBasicCrud(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_chars.DnaTestChecker.build()

    def setUp(self):
        self.expected_taxon_namespaces = []
        self.standalone_taxon_namespaces = []
        self.standalone_taxon_namespaces.append(dendropy.TaxonNamespace(["t1", "t2", "t3"]))
        self.expected_taxon_namespaces.append(self.standalone_taxon_namespaces[0])
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
            pdo1 = standard_file_test_chars.DnaTestChecker.get_char_matrix()
            self.expected_char_matrices[pdo1] = pdo1.taxon_namespace
            self.expected_taxon_namespaces.append(pdo1.taxon_namespace)
            for j in range(2):
                pdo2 = standard_file_test_chars.ProteinTestChecker.get_char_matrix(taxon_namespace=pdo1.taxon_namespace)
                self.expected_char_matrices[pdo2] = pdo2.taxon_namespace

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

        # ds = dendropy.DataSet(label="d1")
        # tree_list =
        # trees = [t for t in tree_list]
        # taxon_namespace = tree_list.taxon_namespace
        # self.assertEqual(len(ds.tree_lists), 0)
        # self.assertEqual(len(ds.taxon_namespaces), 0)
        # ds.add_tree_list(tree_list)
        # self.assertEqual(len(ds.tree_lists), 1)
        # self.assertIs(ds.tree_lists[0], tree_list)
        # self.assertEqual(len(ds.taxon_namespaces), 1)
        # self.assertIs(ds.taxon_namespaces[0], taxon_namespace)
        # self.assertIs(ds.tree_lists[0], taxon_namespace)
        # self.assertEqual(len(ds.tree_lists[0]), len(trees))
        # for t1, t2 in zip(ds.tree_lists[0], trees):
        #     self.assertIs(t1, t2)

if __name__ == "__main__":
    unittest.main()


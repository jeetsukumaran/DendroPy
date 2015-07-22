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

import dendropy
import collections
from dendropy.test.support import curated_test_tree
from dendropy.test.support import standard_file_test_chars

class StandardSingleTaxonNamespaceDataSet(curated_test_tree.CuratedTestTree):

    @staticmethod
    def build(cls):
        standard_file_test_chars.DnaTestChecker.build()
        standard_file_test_chars.RnaTestChecker.build()
        standard_file_test_chars.ProteinTestChecker.build()
        standard_file_test_chars.Standard01234TestChecker.build()

    def verify_chars(self, ds):
        self.assertEqual(len(ds.taxon_namespaces), 1)
        tns = ds.taxon_namespaces[0]
        checkers = (
                standard_file_test_chars.RnaTestChecker,
                standard_file_test_chars.ProteinTestChecker,
                standard_file_test_chars.Standard01234TestChecker,
                standard_file_test_chars.DnaTestChecker,
                )
        self.assertEqual(len(ds.char_matrices), len(checkers))
        for idx, (char_matrix, checker) in enumerate(zip(ds.char_matrices, checkers)):
            self.assertIs(char_matrix.taxon_namespace, tns)
            if checker.matrix_type is dendropy.StandardCharacterMatrix:
                checker.create_class_fixtures_label_sequence_map_based_on_state_alphabet(checker,
                        char_matrix.default_state_alphabet)
            standard_file_test_chars.general_char_matrix_checker(
                    self,
                    char_matrix,
                    checker,
                    check_taxon_annotations=self.check_taxon_annotations,
                    check_matrix_annotations=self.check_matrix_annotations,
                    check_sequence_annotations=self.check_sequence_annotations,
                    check_column_annotations=self.check_column_annotations,
                    check_cell_annotations=self.check_cell_annotations,)

    def verify_trees(self, ds):
        self.assertEqual(len(ds.taxon_namespaces), 1)
        tns = ds.taxon_namespaces[0]
        self.assertEqual(len(ds.tree_lists), 7)
        for tree_list_idx, tree_list in enumerate(ds.tree_lists):
            self.assertIs(tree_list.taxon_namespace, tns)
            expected_labels = (
                    'the first tree',
                    'the SECOND tree',
                    'The Third Tree',
                    )
            self.assertEqual(len(tree_list), len(expected_labels))
            for tree_idx, (tree, expected_label) in enumerate(zip(tree_list, expected_labels)):
                # print(tree_list_idx, tree_idx)
                self.assertEqual(tree.label, expected_label)
                self.verify_curated_tree(
                        tree=tree,
                        suppress_internal_node_taxa=False,
                        suppress_leaf_node_taxa=False,
                        suppress_edge_lengths=False,
                        node_taxon_label_map=None)

    def verify_dataset(self, ds):
        self.verify_chars(ds)
        self.verify_trees(ds)

class MultipleTaxonNamespaceDataSet(object):

    expected_taxa = collections.OrderedDict()
    expected_taxa["X1"] = ("x1.1", "x1.2", "x1.3", "x1.4", "x1.5", )
    expected_taxa["X2"] = ("x2.1", "x2.2", "x2.3", "x2.4", )
    expected_taxa["X3"] = ("x3.1", "x3.2", "x3.3", "x3.4", "x3.5", "x3.6", )
    expected_chars = collections.OrderedDict()
    expected_chars["X1"] = ("x1.chars1", "x1.chars2", "x1.chars3", )
    expected_chars["X2"] = ("x2.chars1", "x2.chars2")
    expected_trees = collections.OrderedDict()
    expected_trees["X1"] = ("x1.trees1", "x1.trees3", "x1.trees3",)
    expected_trees["X2"] = ("x2.trees1", "x2.trees2")

    def verify_attached_taxon_namespace_written(self, ds, taxon_namespace):
        self.assertEqual(len(ds.taxon_namespaces), 1)
        expected_tns_label = taxon_namespace.label
        # self.assertEqual(ds.taxon_namespaces[0].label, expected_tns_label)
        expected_tlsts = MultipleTaxonNamespaceDataSet.expected_trees.get(expected_tns_label, [])
        self.assertEqual(len(ds.char_matrices), len(expected_tlsts))
        expected_cms = MultipleTaxonNamespaceDataSet.expected_chars.get(expected_tns_label, [])
        self.assertEqual(len(ds.tree_lists), len(expected_cms))

    def verify_attached_taxon_namespace(self, ds, attached_taxon_namespace):
        self.assertEqual(len(ds.taxon_namespaces), 1)
        self.assertIs(ds.taxon_namespaces[0], attached_taxon_namespace)
        self.assertEqual(len(ds.taxon_namespaces[0]),
                sum(len(v) for v in MultipleTaxonNamespaceDataSet.expected_taxa.values()))
        self.assertEqual(len(ds.char_matrices),
                sum(len(v) for v in MultipleTaxonNamespaceDataSet.expected_chars.values()))
        self.assertEqual(len(ds.tree_lists),
                sum(len(v) for v in MultipleTaxonNamespaceDataSet.expected_trees.values()))
        for block in (ds.tree_lists, ds.char_matrices):
            for item in block:
                self.assertIs(item.taxon_namespace, attached_taxon_namespace)
                if hasattr(item[0], "taxon_namespace"):
                    for subitem in item:
                        self.assertIs(subitem.taxon_namespace, attached_taxon_namespace)

    def verify_unrestricted(self, ds):
        self.assertEqual(len(ds.taxon_namespaces), len(MultipleTaxonNamespaceDataSet.expected_taxa))
        self.assertEqual(len(ds.char_matrices),
                sum(len(v) for v in MultipleTaxonNamespaceDataSet.expected_chars.values()))
        self.assertEqual(len(ds.tree_lists),
                sum(len(v) for v in MultipleTaxonNamespaceDataSet.expected_trees.values()))
        for tns, expected_tns_label in zip(ds.taxon_namespaces, MultipleTaxonNamespaceDataSet.expected_taxa):
            self.assertEqual(tns.label, expected_tns_label)
            expected_tns = MultipleTaxonNamespaceDataSet.expected_taxa[expected_tns_label]
            self.assertEqual(len(tns), len(expected_tns))
            for taxon, expected_label in zip(tns, expected_tns):
                self.assertEqual(taxon.label, expected_label)
            cms = [char_matrix for char_matrix in ds.char_matrices if char_matrix.taxon_namespace is tns]
            expected_cms = MultipleTaxonNamespaceDataSet.expected_chars.get(expected_tns_label, [])
            self.assertEqual(len(cms), len(expected_cms))
            tlst = [tree_list for tree_list in ds.tree_lists if tree_list.taxon_namespace is tns]
            expected_tlsts = MultipleTaxonNamespaceDataSet.expected_trees.get(expected_tns_label, [])
            self.assertEqual(len(tlst), len(expected_tlsts))




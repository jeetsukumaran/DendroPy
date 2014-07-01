#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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

import dendropy
from dendropy.test.support import curated_test_tree
from dendropy.test.support import standard_file_test_chars

class StandardDataSet(curated_test_tree.CuratedTestTree):

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
                checker.create_class_data_label_sequence_map_based_on_state_alphabet(checker,
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

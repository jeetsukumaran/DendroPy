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
Tests for general NEWICK tree list reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
from dendropy.test import base_newick_test_cases
from dendropy.test import base_standard_trees_parsing_test_cases

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NewickTreeListReaderDefaultTestCase(
        base_standard_trees_parsing_test_cases.StandardTreesParsingTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls)

    # ## NOTE: tests are in base_standard_trees_parsing_test_cases.StandardTreesParsingTestCase, !! ##

class NewickTreeListReaderSuppressedInternalSuppressedLeafTaxaTestCase(
        base_newick_test_cases.NewickTreeListReaderTaxaManagementBaseTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True,
                )

class NewickTreeListReaderUnsuppressedInternalSuppressedLeafTaxaTestCase(
        base_newick_test_cases.NewickTreeListReaderTaxaManagementBaseTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls,
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=True,
                )

class NewickTreeListReaderSuppressedInternalUnsuppressedLeafTaxaTestCase(
        base_newick_test_cases.NewickTreeListReaderTaxaManagementBaseTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False,
                )

class NewickTreeListReaderUnsuppressedInternalUnsuppressedLeafTaxaTestCase(
        base_newick_test_cases.NewickTreeListReaderTaxaManagementBaseTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls,
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False,
                )


class NewickTreeListMetadataTest(
        base_newick_test_cases.NewickTreeListReaderTaxaManagementBaseTestCase,
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls,
                is_metadata_extracted=True,
                )

    def test_read_metadata(self):
        tree_file_titles = [
            "dendropy-test-trees-multifurcating-rooted-annotated",
            "dendropy-test-trees-n33-unrooted-annotated-x10a",
        ]
        for tree_file_title in tree_file_titles:
            tree_filepath = standard_file_test_trees._TREE_FILEPATHS["newick"][tree_file_title]
            with open(tree_filepath, "r") as src:
                tree_string = src.read()
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(src,
                            "newick",
                            extract_comment_metadata=True)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title)

class NewickTreeListRootingAndMetadataTest(
        standard_file_test_trees.StandardTestTreesChecker,
        dendropytest.ExtendedTestCase):

    def test_correct_rooting_weighting_and_metadata_association(self):
        tree_str = """\
                ;;;;
                [&color=red][&W 0.25][&R](a,(b,(c,d)))[&W 0.1][&color=wrong1];
                [&W 0.1][&color=wrong1][&U];[&W 0.1][&color=wrong1];[&W 0.1][&color=wrong1];
                [&color=red][&W 0.25][&R](a,(b,(c,d)))[&W 0.1][&color=wrong1];
                [&W 0.1][&color=wrong1][&U];[&W 0.1][&color=wrong1];[&W 0.1][&color=wrong1];
                (a,(b,(c,d)));;;
        """
        trees = dendropy.TreeList.get_from_string(tree_str,
                "newick",
                extract_comment_metadata=True,
                store_tree_weights=True)
        self.assertEqual(len(trees.taxon_namespace), 4)
        tax_labels = [t.label for t in trees.taxon_namespace]
        self.assertSequenceEqual(set(tax_labels), set(["a", "b", "c", "d"]))
        self.assertEqual(len(trees), 3)
        for tree_idx, tree in enumerate(trees):
            if tree_idx < 2:
                self.assertIs(tree.is_rooted, True)
                self.assertEqual(tree.weight, 0.25)
                self.assertEqual(tree.annotations.get_value("color", None), "red")
            else:
                self.assertIs(tree.is_rooted, None)
                self.assertEqual(tree.weight, 1.0)
                self.assertFalse(tree.has_annotations)

class NewickTreeListReaderMultipleRedundantSemiColons(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def test_multiple_redundant_semicolons(self):
        tree_str = self.get_newick_string()
        s = ";;;;;{tree_str};;; ;\n; \n ;       ;;{tree_str};;;  [(a,(b,c)];  ; ;;".format(tree_str=tree_str)
        trees = dendropy.TreeList.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False,
                suppress_edge_lengths=False)
        self.assertEqual(len(trees), 2)
        for t in trees:
            self.verify_curated_tree(t,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False,
                suppress_edge_lengths=False)

class NewickTreeListReaderTaxonNamespaceTest(dendropytest.ExtendedTestCase):

    def test_shared_taxon_namespace(self):
        tree_filenames = [
                ("pythonidae.reference-trees.newick", 33), # ntax = 33
                ("pythonidae.reference-trees.newick", 33), # ntax = 33
                ("bird_orders.newick", 56), # ntax = 23
                ("pythonidae.reference-trees.taxon-numbers-only.newick", 89), # ntax = 33
                ("pythonidae.reference-trees.newick", 89), # ntax = 33
                ("bird_orders.newick", 89), # ntax = 23
        ]
        common_taxon_namespace = dendropy.TaxonNamespace()
        prev_expected_ntax = 0
        for tree_filename, expected_ntax in tree_filenames:
            self.assertEqual(len(common_taxon_namespace), prev_expected_ntax)
            tree_filepath = pathmap.tree_source_path(tree_filename)
            for reps in range(3):
                tree_list = dendropy.TreeList.get_from_path(
                        pathmap.tree_source_path(tree_filename),
                        "newick",
                        taxon_namespace=common_taxon_namespace)
                self.assertEqual(len(common_taxon_namespace), expected_ntax)
            prev_expected_ntax = expected_ntax

if __name__ == "__main__":
    unittest.main()

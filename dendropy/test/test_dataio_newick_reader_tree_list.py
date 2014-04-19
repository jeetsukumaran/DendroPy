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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import dendropy
import random
import collections
import json
from dendropy.test.support import datagen_standard_file_test_trees
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import pathmap
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NewickTreeListReaderStandardTestTreeTest(
        datagen_standard_file_test_trees.StandardTestTreeChecker,
        unittest.TestCase):

    schema_tree_filepaths = dict(datagen_standard_file_test_trees.tree_filepaths["newick"])

    def test_default_newick_get(self):
        for tree_file_title in [
                'standard-test-trees-n14-unrooted-treeshapes',
                'standard-test-trees-n10-rooted-treeshapes',
                ]:
            tree_filepath = self.schema_tree_filepaths[tree_file_title]
            with open(tree_filepath, "r") as src:
                tree_string = src.read()
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(src, "newick")
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            metadata_extracted=False,
                            distinct_nodes_and_edges=False)

    def test_default_newick_read(self):
        preloaded_tree_file_title = "standard-test-trees-n33-x10a"
        preloaded_tree_reference = datagen_standard_file_test_trees.tree_references[preloaded_tree_file_title]
        tree_file_title = "standard-test-trees-n33-x10a"
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                # prepopulate
                tree_list = dendropy.TreeList.get_from_path(
                        self.schema_tree_filepaths[preloaded_tree_file_title],
                        "newick")
                # check to make sure trees were loaded
                old_len = len(tree_list)
                self.assertEqual(old_len, len(tree_list._trees))
                self.assertEqual(old_len, preloaded_tree_reference["num_trees"])
                self.verify_standard_trees(
                        tree_list,
                        preloaded_tree_file_title,
                        distinct_nodes_and_edges=False)

                # load
                old_id = id(tree_list)
                f = getattr(tree_list, method)
                trees_read = f(src, "newick")
                new_id = id(tree_list)
                self.assertEqual(old_id, new_id)

                # make sure new trees added
                new_len = len(tree_list)
                self.assertEqual(new_len, len(tree_list._trees))
                expected_number_of_trees = tree_reference["num_trees"]
                self.assertEqual(old_len + expected_number_of_trees, new_len)
                self.assertEqual(trees_read, expected_number_of_trees)

                # check new trees
                for tree_idx, tree in enumerate(tree_list[old_len:]):
                    self.compare_to_check_tree(
                            tree=tree,
                            tree_file_title=tree_file_title,
                            check_tree_idx=tree_idx,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            metadata_extracted=False,
                            distinct_nodes_and_edges=False)

                # make sure old ones still intact
                for tree_idx, tree in enumerate(tree_list[:old_len]):
                    self.compare_to_check_tree(
                            tree=tree,
                            tree_file_title=preloaded_tree_file_title,
                            check_tree_idx=tree_idx,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            metadata_extracted=False,
                            distinct_nodes_and_edges=False)

    def test_selective_taxa_newick_get(self):
        # skip big files
        tree_file_title = "standard-test-trees-n12-x2"
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for suppress_internal_node_taxa in [True, False]:
            for suppress_external_node_taxa in [True, False]:
                kwargs = {
                        "suppress_internal_node_taxa": suppress_internal_node_taxa,
                        "suppress_external_node_taxa": suppress_external_node_taxa,
                }
                with open(tree_filepath, "r") as tree_stream:
                    approaches = (
                            (dendropy.TreeList.get_from_path, tree_filepath),
                            (dendropy.TreeList.get_from_stream, tree_stream),
                            (dendropy.TreeList.get_from_string, tree_string),
                            )
                    for method, src in approaches:
                        tree_list = method(src, "newick", **kwargs)
                        self.verify_standard_trees(
                                tree_list=tree_list,
                                tree_file_title=tree_file_title,
                                suppress_internal_node_taxa=suppress_internal_node_taxa,
                                suppress_external_node_taxa=suppress_external_node_taxa,
                                metadata_extracted=False,
                                distinct_nodes_and_edges=False)

    def test_selective_taxa_newick_read(self):
        # skip big files
        tree_file_title = "standard-test-trees-n12-x2"
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for suppress_internal_node_taxa in [True, False]:
            for suppress_external_node_taxa in [True, False]:
                kwargs = {
                        "suppress_internal_node_taxa": suppress_internal_node_taxa,
                        "suppress_external_node_taxa": suppress_external_node_taxa,
                }
                with open(tree_filepath, "r") as tree_stream:
                    approaches = (
                            ("read_from_path", tree_filepath),
                            ("read_from_stream", tree_stream),
                            ("read_from_string", tree_string),
                            )
                    for method, src in approaches:
                        tree_list = dendropy.TreeList()
                        old_id = id(tree_list)
                        f = getattr(tree_list, method)
                        f(src, "newick", **kwargs)
                        new_id = id(tree_list)
                        self.verify_standard_trees(
                                tree_list=tree_list,
                                tree_file_title=tree_file_title,
                                suppress_internal_node_taxa=suppress_internal_node_taxa,
                                suppress_external_node_taxa=suppress_external_node_taxa,
                                metadata_extracted=False,
                                distinct_nodes_and_edges=False)

    def test_tree_offset_newick_get(self):
        tree_file_title = "standard-test-trees-n33-x100a"
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(
                            src,
                            "newick",
                            collection_offset=0,
                            tree_offset=tree_offset,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            distinct_nodes_and_edges=False)

    def test_tree_offset_newick_read(self):
        tree_file_title = "standard-test-trees-n33-x100a"
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        ("read_from_path", tree_filepath),
                        ("read_from_stream", tree_stream),
                        ("read_from_string", tree_string),
                        )
                for method, src in approaches:
                    tree_list = dendropy.TreeList()
                    f = getattr(tree_list, method)
                    trees_read = f(src,
                            "newick",
                            collection_offset=0,
                            tree_offset=tree_offset,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            distinct_nodes_and_edges=False)

    def test_tree_offset_without_collection_offset_newick_get(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        approaches = (
                dendropy.TreeList.get_from_path,
                dendropy.TreeList.get_from_stream,
                dendropy.TreeList.get_from_string,
                )
        for approach in approaches:
            with self.assertRaises(TypeError):
                approach(tree_filepath, "newick", collection_offset=None, tree_offset=0)

    def test_tree_offset_without_collection_offset_newick_read(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        approaches = (
                "read_from_path",
                "read_from_stream",
                "read_from_string",
                )
        for approach in approaches:
            tree_list = dendropy.TreeList()
            f = getattr(tree_list, approach)
            with self.assertRaises(TypeError):
                f(tree_filepath, "newick", collection_offset=None, tree_offset=0)

    def test_out_of_range_tree_offset_newick_get(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, "newick", collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_tree_offset_newick_read(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, "newick", collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_collection_offset_newick_get(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, "newick", collection_offset=1, tree_offset=0)

    def test_out_of_range_collection_offset_newick_read(self):
        tree_file_title = 'standard-test-trees-n33-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, "newick", collection_offset=1, tree_offset=0)

    def test_unsupported_keyword_arguments_newick_get(self):
        tree_file_title = 'standard-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(TypeError):
                    method(src,
                            "newick",
                            suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                            gobbledegook=False,
                            )

    def test_unsupported_keyword_arguments_newick_read(self):
        tree_file_title = 'standard-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(TypeError):
                    f(src,
                      "newick",
                      suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                      gobbledegook=False,
                    )

class NewickTreeListReaderMultipleRedundantSemiColons(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def test_multiple_redundant_semicolons(self):
        tree_str = self.get_newick_string()
        s = ";;;;;{tree_str};;; ;\n; \n ;       ;;{tree_str};;;  [(a,(b,c)];  ; ;;".format(tree_str=tree_str)
        trees = dendropy.TreeList.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=False,
                suppress_edge_lengths=False)
        self.assertEqual(len(trees), 2)
        for t in trees:
            self.verify_curated_tree(t,
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=False,
                suppress_edge_lengths=False)

class NewickTreeListReaderTaxonNamespaceTest(unittest.TestCase):

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

class NewickTreeListMetadataTest(
        datagen_standard_file_test_trees.StandardTestTreeChecker,
        unittest.TestCase):

    def test_read_metadata(self):
        tree_file_titles = [
                'standard-test-trees-n33-annotated',
        ]
        for tree_file_title in tree_file_titles:
            tree_filepath = datagen_standard_file_test_trees.tree_filepaths["newick"][tree_file_title]
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
                            tree_file_title=tree_file_title,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            metadata_extracted=True,
                            distinct_nodes_and_edges=False)


if __name__ == "__main__":
    unittest.main()

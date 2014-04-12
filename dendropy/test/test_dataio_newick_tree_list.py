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
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
_LOG = get_logger(__name__)

class NewickTreeListReaderStandardTestTreeTest(unittest.TestCase):

    schema_tree_filepaths = dict(datagen_standard_file_test_trees.tree_filepaths["newick"])

    def compare_comments(self,
            item,
            check,
            metadata_extracted=False):
        check_comments = list(check["comments"])
        item_comments = list(item.comments)
        for comment in item.comments:
            try:
                check_comments.remove(comment)
            except ValueError:
                pass
            else:
                item_comments.remove(comment)
        self.assertEqual(check_comments, [])
        if metadata_extracted:
            self.assertEqual(item_comments, [])
        else:
            for idx, c in enumerate(item_comments):
                if c.startswith("&"):
                    item_comments[idx] = c[1:]
            item_metadata_comments = ",".join(item_comments)
            check_metadata_comments = ",".join(check["metadata_comments"])
            self.maxDiff = None
            self.assertEqual(item_metadata_comments, check_metadata_comments)

    def label_nodes(self, tree):
        for node_idx, node in enumerate(tree):
            if node.taxon is not None:
                node.canonical_label = node.taxon.label
            else:
                node.canonical_label = node.label

    def compare_to_check_tree(self,
            tree,
            tree_file_title,
            check_tree_idx,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False,
            metadata_extracted=False):
        check_tree = datagen_standard_file_test_trees.tree_references[tree_file_title][str(check_tree_idx)]
        self.assertIs(tree.is_rooted, check_tree["is_rooted"])
        self.compare_comments(tree, check_tree, metadata_extracted)
        seen_taxa = []
        node_labels = []
        edge_labels = []
        num_visited_nodes = 0
        self.label_nodes(tree)
        for node_idx, node in enumerate(tree):
            num_visited_nodes += 1
            check_node = check_tree["nodes"][node.canonical_label]
            check_node_label = check_node["label"]
            self.assertEqual(node.canonical_label, check_node_label)
            # node_labels.append(node.canonical_label)
            _LOG.debug("{}: {}: {}".format(tree_file_title, check_tree_idx, node.canonical_label))

            check_node_children = check_node["children"]
            if check_node_children:
                self.assertTrue(node.is_internal())
                self.assertFalse(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node_children))
                if suppress_internal_node_taxa:
                    self.assertEqual(node.label, check_node_label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node_label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)
            else:
                self.assertFalse(node.is_internal())
                self.assertTrue(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node_children))
                if suppress_external_node_taxa:
                    self.assertEqual(node.label, check_node_label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node_label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)

            if node.parent_node is not None:
                if node.parent_node.is_internal:
                    if suppress_internal_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node["parent"])
                        self.assertIs(node.parent_node.label, None)
                else:
                    if suppress_external_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node["parent"])
                        self.assertIs(node.parent_node.label, None)
            else:
                self.assertEqual(check_node["parent"], "None")

            child_labels = []
            for ch in node.child_node_iter():
                if ch.is_internal():
                    if suppress_internal_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        self.assertIsNot(ch.taxon, None)
                        child_labels.append(ch.taxon.label)
                        self.assertIs(ch.label, None)
                else:
                    if suppress_external_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        self.assertIsNot(ch.taxon, None)
                        child_labels.append(ch.taxon.label)
                        self.assertIs(ch.label, None)
            self.assertEqual(len(child_labels), len(check_node["children"]))
            self.assertEqual(set(child_labels), set(check_node["children"]))

            edge = node.edge
            check_edge = check_tree["edges"][node.canonical_label]
            if edge.tail_node is None:
                self.assertEqual(check_edge["tail_node"], "None")
            else:
                self.assertEqual(edge.tail_node.canonical_label, check_edge["tail_node"])
            self.assertEqual(edge.head_node.canonical_label, check_edge["head_node"])
            self.assertAlmostEqual(node.edge.length, float(check_edge["length"]))

            # This hackery because NEWICK/NEXUS cannot distinguish between
            # node and edge comments, and everything gets lumped in as a
            # node comment
            node.comments += edge.comments
            d = {
                    "comments": check_node["comments"] + check_edge["comments"],
                    "metadata_comments": check_node["metadata_comments"] + check_edge["metadata_comments"],
                    }
            self.compare_comments(node, d, metadata_extracted)

        self.assertEqual(num_visited_nodes, len(check_tree["nodeset"]))
        self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
        self.assertEqual(set(seen_taxa), set(tree.taxon_namespace))
        node_labels.extend([t.label for t in tree.taxon_namespace])
        self.assertEqual(len(node_labels), len(check_tree["nodeset"]))
        self.assertEqual(set(node_labels), set(check_tree["nodeset"]))

    def verify_standard_trees(self,
            tree_list,
            tree_file_title,
            tree_offset=0,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False,
            metadata_extracted=False):
        tree_reference = datagen_standard_file_test_trees.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        if tree_offset < 0:
            if abs(tree_offset) > expected_number_of_trees:
                tree_offset = 0
            else:
                tree_offset = expected_number_of_trees + tree_offset
        self.assertEqual(len(tree_list), expected_number_of_trees-tree_offset)
        # for tree_idx, (tree, check_tree) in enumerate(zip(tree_list, datagen_standard_file_test_trees.tree_directory[tree_file_title])):
        for tree_idx, tree in enumerate(tree_list):
            _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            self.compare_to_check_tree(
                    tree=tree,
                    tree_file_title=tree_file_title,
                    check_tree_idx=tree_idx + tree_offset,
                    suppress_internal_node_taxa=suppress_internal_node_taxa,
                    suppress_external_node_taxa=suppress_external_node_taxa,
                    metadata_extracted=metadata_extracted)

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
                            metadata_extracted=False)

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
                self.verify_standard_trees(tree_list, preloaded_tree_file_title)

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
                            metadata_extracted=False)

                # make sure old ones still intact
                for tree_idx, tree in enumerate(tree_list[:old_len]):
                    self.compare_to_check_tree(
                            tree=tree,
                            tree_file_title=preloaded_tree_file_title,
                            check_tree_idx=tree_idx,
                            suppress_internal_node_taxa=True,
                            suppress_external_node_taxa=False,
                            metadata_extracted=False)

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
                                metadata_extracted=False)

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
                                metadata_extracted=False)

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
                            suppress_external_node_taxa=False)

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
                            suppress_external_node_taxa=False)

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

    # def test_out_of_range_collection_offset_newick_get_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     t = dendropy.TreeList()
    #     for collection_offset in [1, -2]:
    #         with self.assertRaises(IndexError):
    #             t = dendropy.TreeList.get_from_path(
    #                     pathmap.tree_source_path(tree_filename),
    #                     "newick",
    #                     collection_offset=1,
    #                     tree_offset=0)

    # def test_out_of_range_collection_offset_newick_read_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     t = dendropy.TreeList()
    #     for collection_offset in [1, -2]:
    #         with self.assertRaises(IndexError):
    #             t.read_from_path(
    #                     pathmap.tree_source_path(tree_filename),
    #                     "newick",
    #                     collection_offset=1,
    #                     tree_offset=0)

    # def test_invalid_tree_offset_newick_get_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     with self.assertRaises(IndexError):
    #         t = dendropy.TreeList.get_from_path(
    #                 pathmap.tree_source_path(tree_filename),
    #                 "newick",
    #                 collection_offset=0,
    #                 tree_offset=-1)

    # def test_invalid_tree_offset_newick_read_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     t = dendropy.TreeList()
    #     with self.assertRaises(IndexError):
    #         t.read_from_path(
    #                 pathmap.tree_source_path(tree_filename),
    #                 "newick",
    #                 collection_offset=0,
    #                 tree_offset=-1)

    # def test_invalid_collection_offset_newick_get_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     t = dendropy.TreeList()
    #     with self.assertRaises(IndexError):
    #         t = dendropy.TreeList.get_from_path(
    #                 pathmap.tree_source_path(tree_filename),
    #                 "newick",
    #                 collection_offset=-1,
    #                 tree_offset=0)

    # def test_invalid_collection_offset_newick_read_from(self):
    #     tree_filename = datagen_standard_file_test_trees.newick_tree_filenames[0]
    #     tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
    #     t = dendropy.TreeList()
    #     with self.assertRaises(IndexError):
    #         n = t.read_from_path(
    #                 pathmap.tree_source_path(tree_filename),
    #                 "newick",
    #                 collection_offset=-1,
    #                 tree_offset=0)

# class NewickTreeListReaderTaxonNamespaceTest(unittest.TestCase):

#     def test_shared_taxon_namespace(self):
#         tree_filenames = [
#             ("pythonidae.reference-trees.newick", 33), # ntax = 33
#             ("pythonidae.reference-trees.newick", 33), # ntax = 33
#             ("bird_orders.newick", 56), # ntax = 23
#             ("pythonidae.reference-trees.taxon-numbers-only.newick", 89), # ntax = 33
#             ("pythonidae.reference-trees.newick", 89), # ntax = 33
#             ("bird_orders.newick", 89), # ntax = 23
#             ]
#         common_taxon_namespace = dendropy.TaxonNamespace()
#         prev_expected_ntax = 0
#         for tree_filename, expected_ntax in tree_filenames:
#             self.assertEqual(len(common_taxon_namespace), prev_expected_ntax)
#             tree_filepath = pathmap.tree_source_path(tree_filename)
#             for reps in range(3):
#                 tree_list = dendropy.TreeList.get_from_path(
#                         pathmap.tree_source_path(tree_filename),
#                         "newick",
#                         taxon_namespace=common_taxon_namespace)
#                 self.assertEqual(len(common_taxon_namespace), expected_ntax)
#             prev_expected_ntax = expected_ntax

if __name__ == "__main__":
    unittest.main()

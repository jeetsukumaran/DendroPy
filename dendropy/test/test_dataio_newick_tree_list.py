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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import dendropy
import random
from dendropy.test.support import standard_test_tree_data
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
_LOG = get_logger(__name__)

class NewickTreeListReaderStandardTestTreeTest(unittest.TestCase):

    def compare_to_check_tree(self,
            tree,
            tree_file_title,
            check_tree_idx,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False):
        check_tree = standard_test_tree_data.tree_directory[tree_file_title][check_tree_idx]
        self.assertEqual(tree.comments, check_tree.comments)
        self.assertIs(tree.is_rooted, check_tree.rooted)
        seen_taxa = []
        node_labels = []
        num_visited_nodes = 0
        for node_idx, node in enumerate(tree):
            num_visited_nodes += 1
            if node.taxon is not None:
                check_node_label_key = node.taxon.label
            else:
                check_node_label_key = node.label
            check_node = standard_test_tree_data.node_directory[(tree_file_title, check_tree_idx, check_node_label_key)]
            _LOG.debug("{}: {}: {}".format(tree_file_title, check_tree_idx, check_node_label_key))
            if check_node.children:
                self.assertTrue(node.is_internal())
                self.assertFalse(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node.children))
                if suppress_internal_node_taxa:
                    self.assertEqual(node.label, check_node.label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node.label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)
            else:
                self.assertFalse(node.is_internal())
                self.assertTrue(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node.children))
                if suppress_external_node_taxa:
                    self.assertEqual(node.label, check_node.label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node.label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)
            self.assertAlmostEqual(node.edge.length, check_node.edge_length, 4)
            if node.parent_node is not None:
                if node.parent_node.is_internal:
                    if suppress_internal_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node.parent)
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node.parent)
                        self.assertIs(node.parent_node.label, None)
                else:
                    if suppress_external_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node.parent)
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node.parent)
                        self.assertIs(node.parent_node.label, None)
            else:
                self.assertEqual(check_node.parent, "None")
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
            self.assertEqual(len(child_labels), len(check_node.children))
            self.assertEqual(set(child_labels), set(check_node.children))
            self.assertEqual(node.comments, check_node.comments)
        self.assertEqual(num_visited_nodes, len(check_tree.nodes))
        self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
        self.assertEqual(set(seen_taxa), set(tree.taxon_namespace))
        node_labels.extend([t.label for t in tree.taxon_namespace])
        self.assertEqual(len(node_labels), len(check_tree.nodes))
        self.assertEqual(set(node_labels), set(check_tree.nodes))

    def verify_standard_trees(self,
            tree_list,
            tree_file_title,
            tree_offset=0,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False):
        self.assertEqual(len(tree_list), standard_test_tree_data.expected_number_of_trees[tree_file_title]-tree_offset)
        # for tree_idx, (tree, check_tree) in enumerate(zip(tree_list, standard_test_tree_data.tree_directory[tree_file_title])):
        for tree_idx, tree in enumerate(tree_list):
            _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            self.compare_to_check_tree(
                    tree=tree,
                    tree_file_title=tree_file_title,
                    check_tree_idx=tree_idx + tree_offset,
                    suppress_internal_node_taxa=suppress_internal_node_taxa,
                    suppress_external_node_taxa=suppress_external_node_taxa)

    def test_default_newick_get_from_path(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick")
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_default_newick_get_from_string(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            filepath = pathmap.tree_source_path(tree_filename)
            fsrc = open(filepath, "r", newline=None)
            with fsrc:
                s = fsrc.read()
            tree_list = dendropy.TreeList.get_from_string(s, "newick")
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_default_newick_get_from_stream(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            filepath = pathmap.tree_source_path(tree_filename)
            fsrc = open(filepath, "r", newline=None)
            with fsrc:
                tree_list = dendropy.TreeList.get_from_stream(fsrc, "newick")
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_default_newick_read_from_path(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            n = tree_list.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick")
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_default_newick_read_from_string(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            filepath = pathmap.tree_source_path(tree_filename)
            fsrc = open(filepath, "r", newline=None)
            with fsrc:
                s = fsrc.read()
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            n = tree_list.read_from_string(s, "newick")
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_default_newick_read_from_stream(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            filepath = pathmap.tree_source_path(tree_filename)
            fsrc = open(filepath, "r", newline=None)
            with fsrc:
                n = tree_list.read_from_stream(fsrc, "newick")
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_notaxa_newick_get_from(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList()
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)

    def test_notaxa_newick_read_from(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            n = tree_list.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)

    def test_all_taxa_newick_get_from(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False)
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False)

    def test_all_taxa_newick_read_from(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            n = tree_list.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False)
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    suppress_internal_node_taxa=False,
                    suppress_external_node_taxa=False)

    def test_tree_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        offsets = [0, standard_test_tree_data.expected_number_of_trees[tree_file_title]-1]
        for rep in range(5):
            offsets.append(random.randint(1, standard_test_tree_data.expected_number_of_trees[tree_file_title]-1))
        for tree_offset in offsets:
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
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

    def test_tree_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        offsets = [0, standard_test_tree_data.expected_number_of_trees[tree_file_title]-1]
        for rep in range(5):
            offsets.append(random.randint(1, standard_test_tree_data.expected_number_of_trees[tree_file_title]-1))
        for tree_offset in offsets:
            tree_list = dendropy.TreeList()
            old_id = id(tree_list)
            n = tree_list.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=0,
                    tree_offset=tree_offset,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)
            self.assertEqual(id(tree_list), old_id)
            self.assertEqual(n, len(tree_list))
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    tree_offset=tree_offset,
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=False)

    def test_tree_offset_without_collection_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        with self.assertRaises(TypeError):
            t = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=None,
                    tree_offset=0)

    def test_tree_offset_without_collection_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(TypeError):
            t.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=None,
                    tree_offset=0)

    def test_out_of_range_tree_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        with self.assertRaises(IndexError):
            t = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=0,
                    tree_offset=standard_test_tree_data.expected_number_of_trees[tree_file_title])

    def test_out_of_range_tree_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            t.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=0,
                    tree_offset=standard_test_tree_data.expected_number_of_trees[tree_file_title])

    def test_out_of_range_collection_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            t = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=1,
                    tree_offset=0)

    def test_out_of_range_collection_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            t.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=1,
                    tree_offset=0)

    def test_invalid_tree_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        with self.assertRaises(IndexError):
            t = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=0,
                    tree_offset=-1)

    def test_invalid_tree_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            t.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=0,
                    tree_offset=-1)

    def test_invalid_collection_offset_newick_get_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            t = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=-1,
                    tree_offset=0)

    def test_invalid_collection_offset_newick_read_from(self):
        tree_filename = standard_test_tree_data.newick_tree_filenames[0]
        tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
        t = dendropy.TreeList()
        with self.assertRaises(IndexError):
            n = t.read_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    collection_offset=-1,
                    tree_offset=0)

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

if __name__ == "__main__":
    unittest.main()

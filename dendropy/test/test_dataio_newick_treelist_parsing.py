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
from dendropy.test.support import standard_test_tree_data
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class NewickTreeListReaderStandardTestTreeTest(unittest.TestCase):


    def test_normal_newick_reader(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick")
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.assertEqual(len(tree_list), standard_test_tree_data.expected_number_of_trees[tree_file_title])
            for tree_idx, (tree, check_tree) in enumerate(zip(tree_list, standard_test_tree_data.tree_directory[tree_file_title])):
                _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
                self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
                self.assertEqual(tree.comments, check_tree.comments)
                self.assertIs(tree.is_rooted, check_tree.rooted)
                seen_taxa = []
                node_labels = []
                for node_idx, node in enumerate(tree):
                    _LOG.debug("{}: {}: {}".format(tree_file_title, tree_idx, node.label))
                    if node.taxon is not None:
                        check_node_label_key = node.taxon.label
                    else:
                        check_node_label_key = node.label
                    check_node = standard_test_tree_data.node_directory[(tree_file_title, tree_idx, check_node_label_key)]
                    if check_node.children:
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
                        self.assertEqual(node.parent_node.label, check_node.parent)
                    else:
                        self.assertEqual(check_node.parent, "None")
                    child_labels = []
                    for ch in node.child_node_iter():
                        if not ch.is_leaf():
                            self.assertIs(ch.taxon, None)
                            child_labels.append(ch.label)
                        else:
                            self.assertIsNot(ch.taxon, None)
                            child_labels.append(ch.taxon.label)
                            self.assertIs(ch.label, None)
                    self.assertEqual(len(child_labels), len(check_node.children))
                    self.assertEqual(set(child_labels), set(check_node.children))
                    self.assertEqual(node.comments, check_node.comments)
                self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
                self.assertEqual(set(seen_taxa), set(tree.taxon_namespace))
                node_labels.extend([t.label for t in tree.taxon_namespace])
                self.assertEqual(len(node_labels), len(check_tree.nodes))
                self.assertEqual(set(node_labels), set(check_tree.nodes))

    def test_notaxa_newick_reader(self):
        taxon_namespaces = set()
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            tree_list = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)
            self.assertNotIn(tree_list.taxon_namespace, taxon_namespaces)
            taxon_namespaces.add(tree_list.taxon_namespace)
            self.assertEqual(len(tree_list), standard_test_tree_data.expected_number_of_trees[tree_file_title])
            for tree_idx, (tree, check_tree) in enumerate(zip(tree_list, standard_test_tree_data.tree_directory[tree_file_title])):
                _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
                self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
                self.assertEqual(tree.comments, check_tree.comments)
                self.assertIs(tree.is_rooted, check_tree.rooted)
                node_labels = []
                for node_idx, node in enumerate(tree):
                    _LOG.debug("{}: {}: {}".format(tree_file_title, tree_idx, node.label))
                    check_node = standard_test_tree_data.node_directory[(tree_file_title, tree_idx, node.label)]
                    assert check_node.label == node.label
                    self.assertAlmostEqual(node.edge.length, check_node.edge_length, 4)
                    if node.parent_node is not None:
                        self.assertEqual(node.parent_node.label, check_node.parent)
                    else:
                        self.assertEqual(check_node.parent, "None")
                    ch_labels = [ch.label for ch in node.child_node_iter()]
                    self.assertEqual(len(ch_labels), len(check_node.children))
                    self.assertEqual(set(ch_labels), set(check_node.children))
                    self.assertEqual(node.comments, check_node.comments)
                    node_labels.append(node.label)
                self.assertEqual(len(node_labels), len(check_tree.nodes))
                self.assertEqual(set(node_labels), set(check_tree.nodes))
                self.assertEqual(len(tree_list.taxon_namespace), 0)

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

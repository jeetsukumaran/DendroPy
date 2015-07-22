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
Tests for NEXUS tree writing.
"""

import collections
import unittest
import dendropy
import re
from dendropy.test.support import pathmap
from dendropy.test.support import compare_and_validate
from dendropy.test.support import dendropytest
from dendropy.test.support import curated_test_tree
from dendropy.test.support import standard_file_test_trees
from dendropy.test.test_dataio_newick_writer import newick_tree_writer_test_tree

class NexusTreeWriterTests(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    def test_simple(self):
        tree1, all_nodes, leaf_nodes, internal_nodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False
                )
        kwargs = {
                "suppress_leaf_node_labels": True,
                "suppress_internal_node_labels": True
                }
        s = self.write_out_validate_equal_and_return(
                tree1, "nexus", kwargs)
        tree2 = dendropy.Tree.get_from_string(
                s, "nexus",
                )
        self.verify_curated_tree(tree2)

class NexusTreeWriterDefaultTest(
        standard_file_test_trees.NexusTestTreesChecker,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexusTestTreesChecker.create_class_fixtures(cls,
                is_metadata_extracted=True,
                )

    def test_roundtrip_full(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-annotated-x10a'
        tree_filepath = standard_file_test_trees._TREE_FILEPATHS["nexus"][tree_file_title]
        tree1 = dendropy.Tree.get_from_path(
                tree_filepath,
                "nexus",
                extract_comment_metadata=True,
                store_tree_weights=True,
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False,
        )
        kwargs = {
            "suppress_leaf_taxon_labels"     :  False , # default: False ,
            "suppress_leaf_node_labels"      :  True  , # default: True  ,
            "suppress_internal_taxon_labels" :  False , # default: False ,
            "suppress_internal_node_labels"  :  True  , # default: False ,
            "suppress_rooting"               :  False , # default: False ,
            "suppress_edge_lengths"          :  False , # default: False ,
            "unquoted_underscores"           :  False , # default: False ,
            "preserve_spaces"                :  False , # default: False ,
            "store_tree_weights"             :  False , # default: False ,
            "suppress_annotations"           :  False , # default: True  ,
            "annotations_as_nhx"             :  False , # default: False ,
            "suppress_item_comments"         :  False , # default: True  ,
            "node_label_element_separator"   :  ' '   , # default: ' '   ,
            "node_label_compose_fn"        :  None  , # default: None  ,
            "edge_label_compose_fn"        :  None  , # default: None  ,
        }
        s = self.write_out_validate_equal_and_return(
                tree1, "nexus", kwargs)
        tree2 = dendropy.Tree.get_from_string(
                s,
                "nexus",
                extract_comment_metadata=True,
                store_tree_weights=True,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False,
        )
        self.compare_to_reference_by_title_and_index(
            tree=tree2,
            tree_file_title=tree_file_title,
            reference_tree_idx=0)
            # suppress_internal_node_taxa=False,
            # suppress_leaf_node_taxa=False,
            # is_metadata_extracted=True,
            # is_coerce_metadata_values_to_string=True,
            # is_distinct_nodes_and_edges_representation=False)

class NexusTreeWriterGeneralOptionsTests(
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    def test_node_labeling(self):
        for has_leaf_node_taxa in (True, False):
            for has_leaf_node_labels in (True, False):
                for has_internal_node_taxa in (True, False):
                    for has_internal_node_labels in (True, False):
                        for label_separator in (' ', '$$$'):
                            tree = newick_tree_writer_test_tree(
                                    has_leaf_node_taxa=has_leaf_node_taxa,
                                    has_leaf_node_labels=has_leaf_node_labels,
                                    has_internal_node_taxa=has_internal_node_taxa,
                                    has_internal_node_labels=has_internal_node_labels,
                                    label_separator=label_separator,
                                    )
                            for suppress_leaf_taxon_labels in (True, False):
                                for suppress_leaf_node_labels in (True, False):
                                    for suppress_internal_taxon_labels in (True, False):
                                        for suppress_internal_node_labels in (True, False):
                                            kwargs = {
                                                    "suppress_leaf_taxon_labels"     : suppress_leaf_taxon_labels,
                                                    "suppress_leaf_node_labels"      : suppress_leaf_node_labels,
                                                    "suppress_internal_taxon_labels" : suppress_internal_taxon_labels,
                                                    "suppress_internal_node_labels"  : suppress_internal_node_labels,
                                                    "node_label_element_separator"   : label_separator,
                                                    }
                                            s = self.write_out_validate_equal_and_return(
                                                    tree, "nexus", kwargs)
                                            tree2 = dendropy.Tree.get_from_string(
                                                    s,
                                                    "nexus",
                                                    extract_comment_metadata=True,
                                                    store_tree_weights=True,
                                                    suppress_internal_node_taxa=True,
                                                    suppress_leaf_node_taxa=True,
                                            )
                                            nodes1 = [nd for nd in tree]
                                            nodes2 = [nd for nd in tree2]
                                            self.assertEqual(len(nodes1), len(nodes2))
                                            for nd1, nd2 in zip(nodes1, nodes2):
                                                is_leaf = nd1.is_leaf()
                                                self.assertEqual(nd2.is_leaf(), is_leaf)
                                                if is_leaf:
                                                    self.assertEqual(nd2.label,
                                                            nd1.expected_label[(suppress_leaf_taxon_labels, suppress_leaf_node_labels)])
                                                else:
                                                    self.assertEqual(nd2.label,
                                                            nd1.expected_label[ (suppress_internal_taxon_labels, suppress_internal_node_labels) ])

    def test_rooting_token(self):
        tree1 = newick_tree_writer_test_tree()
        for rooted_state in (None, True, False):
            tree1.is_rooted = rooted_state
            for suppress_rooting in (True, False):
                kwargs = {
                        "suppress_rooting": suppress_rooting,
                }
                s = self.write_out_validate_equal_and_return(
                        tree1, "nexus", kwargs)
                tree2 = dendropy.Tree.get_from_string(
                        s, "nexus", rooting=None)
                if suppress_rooting:
                    self.assertTrue(tree2.is_rootedness_undefined)
                else:
                    if rooted_state is True:
                        self.assertTrue(tree2.is_rooted)
                        self.assertFalse(tree2.is_unrooted)
                    elif rooted_state is False:
                        self.assertFalse(tree2.is_rooted)
                        self.assertTrue(tree2.is_unrooted)
                    else:
                        self.assertTrue(tree2.is_rootedness_undefined)

    def test_edge_lengths(self):
        tree1 = newick_tree_writer_test_tree()
        for suppress_edge_lengths in (True, False):
            kwargs = {
                    "suppress_edge_lengths": suppress_edge_lengths,
            }
            s = self.write_out_validate_equal_and_return(
                    tree1, "nexus", kwargs)
            tree2 = dendropy.Tree.get_from_string(
                    s, "nexus", rooting=None)
            nodes1 = [nd for nd in tree1]
            nodes2 = [nd for nd in tree2]
            self.assertEqual(len(nodes1), len(nodes2))
            for nd1, nd2 in zip(nodes1, nodes2):
                if suppress_edge_lengths:
                    self.assertIs(nd2.edge.length, None)
                else:
                    self.assertEqual(nd2.edge.length, nd1.edge.length)

    def test_unquoted_underscores(self):
        tree1 = newick_tree_writer_test_tree(
                has_leaf_node_labels=False,
                has_internal_node_labels=False)
        for taxon in tree1.taxon_namespace:
            taxon.label = "{label}_{label}".format(label=taxon.label)
        for unquoted_underscores in (True, False):
            kwargs = {
                    "unquoted_underscores": unquoted_underscores,
            }
            s = self.write_out_validate_equal_and_return(
                    tree1, "nexus", kwargs)
            for preserve_underscores in (True, False):
                tree2 = dendropy.Tree.get_from_string(
                        s,
                        "nexus",
                        suppress_internal_node_taxa=False,
                        preserve_underscores=preserve_underscores)
                nodes1 = [nd for nd in tree1]
                nodes2 = [nd for nd in tree2]
                self.assertEqual(len(nodes1), len(nodes2))
                for nd1, nd2 in zip(nodes1, nodes2):
                    original_label = nd1.taxon.label
                    if unquoted_underscores:
                        if preserve_underscores:
                            expected_label = original_label
                        else:
                            expected_label = original_label.replace("_", " ")
                    else:
                        expected_label = original_label
                    self.assertEqual(nd2.taxon.label, expected_label)

    def test_preserve_spaces(self):
        tree1 = newick_tree_writer_test_tree(
                has_leaf_node_labels=False,
                has_internal_node_labels=False)
        for taxon in tree1.taxon_namespace:
            taxon.label = "{label} {label}".format(label=taxon.label)
        for preserve_spaces in (True, False):
            kwargs = {
                    "preserve_spaces": preserve_spaces,
            }
            s = self.write_out_validate_equal_and_return(
                    tree1, "nexus", kwargs)
            tree2 = dendropy.Tree.get_from_string(
                    s,
                    "nexus",
                    suppress_internal_node_taxa=False,
                    preserve_underscores=True)
            nodes1 = [nd for nd in tree1]
            nodes2 = [nd for nd in tree2]
            self.assertEqual(len(nodes1), len(nodes2))
            for nd1, nd2 in zip(nodes1, nodes2):
                original_label = nd1.taxon.label
                if preserve_spaces:
                    expected_label = original_label
                else:
                    expected_label = original_label.replace(" ", "_")
                self.assertEqual(nd2.taxon.label, expected_label)

    def test_store_tree_weights(self):
        tree1 = newick_tree_writer_test_tree(
                has_leaf_node_labels=False,
                has_internal_node_labels=False)
        for store_tree_weights in (True, False):
            for weight in (None, "23.0", "1/2", 1.0):
                tree1.weight = weight
                kwargs = {
                        "store_tree_weights": store_tree_weights,
                }
                s = self.write_out_validate_equal_and_return(
                        tree1, "nexus", kwargs)
                tree2 = dendropy.Tree.get_from_string(
                        s,
                        "nexus",
                        store_tree_weights=True)
                if store_tree_weights and weight is not None:
                    self.assertTrue("[&W " in s)
                    try:
                        w = float(weight)
                    except ValueError:
                        w = eval("/".join(str(float(w)) for w in weight.split("/")))
                    self.assertEqual(tree2.weight, w)
                else:
                    self.assertFalse("[&W " in s)
                    self.assertEqual(tree2.weight, 1.0) # default weight

    def test_suppress_annotations(self):
        tree1 = dendropy.Tree()
        a1 = tree1.seed_node.new_child()
        a2 = tree1.seed_node.new_child()
        tree1.annotations.add_new("t", 1)
        for nd in tree1:
            nd.annotations.add_new("a", 1)
            nd.edge.annotations.add_new("b", 2)
        for suppress_annotations in (True, False):
            for annotations_as_nhx in (True, False):
                kwargs = {
                        "suppress_annotations"   :  suppress_annotations,
                        "annotations_as_nhx"     :  annotations_as_nhx,
                }
                s = self.write_out_validate_equal_and_return(
                        tree1, "nexus", kwargs)
                tree2 = dendropy.Tree.get_from_string(
                        s,
                        "nexus",
                        extract_comment_metadata=True)
                if suppress_annotations:
                    self.assertFalse(tree2.has_annotations)
                    for nd in tree2:
                        self.assertFalse(nd.has_annotations)
                else:
                    if annotations_as_nhx:
                        self.assertEqual(s.count("[&&NHX"), 7)
                        # self.assertEqual(len(re.findall(r"\[&&NHX", s)), 7)
                    else:
                        self.assertEqual(s.count("[&&NHX"), 0)
                        self.assertEqual(s.count("[&"), 7)
                        # self.assertEqual(len(re.findall(r"\[&&NHX", s)), 0)
                        # self.assertEqual(len(re.findall(r"\[&.*?\]", s)), 7)
                    self.assertTrue(tree2.has_annotations)
                    self.assertEqual(tree2.annotations.get_value("t"), '1')
                    for nd in tree2:
                        self.assertTrue(nd.has_annotations)
                        self.assertEqual(nd.annotations.get_value("a"), '1')
                        self.assertEqual(nd.annotations.get_value("b"), '2')

    def test_suppress_item_comments(self):
        tree1 = dendropy.Tree()
        a1 = tree1.seed_node.new_child()
        a2 = tree1.seed_node.new_child()
        tree1.comments.append("t1")
        for nd in tree1:
            nd.comments.append("n1")
            nd.edge.comments.append("e1")
        for suppress_item_comments in (True, False):
            kwargs = {
                    "suppress_item_comments"   :  suppress_item_comments,
            }
            s = self.write_out_validate_equal_and_return(
                    tree1, "nexus", kwargs)
            tree2 = dendropy.Tree.get_from_string(
                    s,
                    "nexus",
                    extract_comment_metadata=False)
            if suppress_item_comments:
                self.assertEqual(tree2.comments, [])
                for nd in tree2:
                    self.assertEqual(nd.comments, [])
                    self.assertEqual(nd.edge.comments, [])
            else:
                self.assertEqual(tree2.comments, ["t1"])
                for nd in tree2:
                    self.assertEqual(nd.comments, ["n1", "e1"])

    def test_node_label_compose_fn(self):
        tree1 = dendropy.Tree()
        a1 = tree1.seed_node.new_child(label="a1")
        a1.taxon = tree1.taxon_namespace.require_taxon("hula")
        a2 = tree1.seed_node.new_child(label="a1")
        a2.taxon = tree1.taxon_namespace.require_taxon("hoop")
        f = lambda x: "zzz"
        kwargs = {
                "suppress_leaf_taxon_labels"     :  False ,
                "suppress_leaf_node_labels"      :  False ,
                "suppress_internal_taxon_labels" :  False ,
                "suppress_internal_node_labels"  :  False ,
                "node_label_compose_fn"   :  f,
        }
        s = self.write_out_validate_equal_and_return(
                tree1, "nexus", kwargs)
        tree2 = dendropy.Tree.get_from_string(
                s,
                "nexus",
                suppress_leaf_node_taxa=True,
                suppress_internal_node_taxa=True)
        for nd in tree2:
            self.assertEqual(nd.label, "zzz")

    def test_edge_label_compose_fn(self):
        tree1 = dendropy.Tree()
        tree1.seed_node.edge.length = 1
        a1 = tree1.seed_node.new_child(label="a1", edge_length=1)
        a2 = tree1.seed_node.new_child(label="a1", edge_length=1)
        f = lambda x: 1000
        kwargs = {
                "edge_label_compose_fn"   :  f,
        }
        s = self.write_out_validate_equal_and_return(
                tree1, "nexus", kwargs)
        tree2 = dendropy.Tree.get_from_string(
                s,
                "nexus",
                suppress_leaf_node_taxa=True,
                suppress_internal_node_taxa=True)
        for nd in tree2:
            self.assertEqual(nd.edge.length, 1000)

if __name__ == "__main__":
    unittest.main()

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
Tests basic Tree structure and iteration.
"""

import unittest
import dendropy
from dendropy.test.support import curated_test_tree
from dendropy.test.support import dendropytest

class TestTreeNodeAndEdgeCollections(curated_test_tree.CuratedTestTree, unittest.TestCase):

    def test_get_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.nodes()
        self.assertEqual(len(nodes), len(anodes))
        self.assertEqual(set(nodes), anodes)
        obs_labels = [nd.label for nd in nodes]

    def test_get_nodes_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.nodes(filter_fn = lambda x : x.edge.length > 13)
        exp_nodes = set([nd for nd in anodes if nd.edge.length > 13])
        for nd in nodes:
            self.assertTrue(nd.edge.length > 13)
        self.assertEqual(len(nodes), len(exp_nodes))
        self.assertEqual(set(nodes), exp_nodes)

    def test_get_leaf_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.leaf_nodes()
        self.assertEqual(len(nodes), len(lnodes))
        self.assertEqual(set(nodes), lnodes)

    def test_get_internal_nodes_with_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.internal_nodes()
        nlnodes = inodes | set([tree.seed_node])
        self.assertEqual(len(nodes), len(nlnodes))
        self.assertEqual(set(nodes), nlnodes)

    def test_get_internal_nodes_no_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.internal_nodes(True)
        self.assertEqual(len(nodes), len(inodes))
        self.assertEqual(set(nodes), inodes)

    def test_get_edges(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.edges()
        eset = set([nd.edge for nd in anodes])
        self.assertEqual(len(edges), len(eset))
        self.assertEqual(set(edges), eset)

    def test_get_edges_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.edges(filter_fn=lambda x : x.length > 13)
        exp_edges = set([nd.edge for nd in anodes if nd.edge.length > 13])
        for edge in edges:
            self.assertTrue(edge.length > 13)
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_leaf_edges(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.leaf_edges()
        exp_edges = set([nd.edge for nd in lnodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_internal_edges_with_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.internal_edges()
        nlnodes = inodes | set([tree.seed_node])
        exp_edges = set([nd.edge for nd in nlnodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_internal_edges_no_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.internal_edges(True)
        exp_edges = set([nd.edge for nd in inodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_child_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for node in anodes:
            child_labels = [ch.label for ch in node.child_nodes()]
            expected_children = self.node_children[node.label]
            self.assertEqual(len(child_labels), len(expected_children))
            self.assertEqual(set(child_labels), set(expected_children))

class TestTreeNodeFinders(curated_test_tree.CuratedTestTree, unittest.TestCase):

    def test_find_node(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node(lambda x: x.label == "c")
        self.assertEqual(node.label, "c")

    def test_find_node_nonexisting(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node(lambda x: x.label == "zzz")
        self.assertIs(node, None)

    def test_find_node_with_label(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node_with_label("c")
        self.assertEqual(node.label, "c")

    def test_find_node_with_label_nonexisting(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node_with_label("zzz")
        self.assertIs(node, None)

class TestTreeIterators(curated_test_tree.CuratedTestTree, unittest.TestCase):

    ### Default Iterator ###

    def test_default_iteration(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.preorder_sequence)

    ### Preorder Node Iterator ###

    def test_preorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.preorder_sequence)

    def test_preorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Preorder Internal Node Iterator ###

    def test_preorder_internal_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_internal_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x]]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_internal_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_internal_node_iter(exclude_seed_node=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_internal_node_iter(exclude_seed_node=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Postorder Node Iterator ###

    def test_postorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.postorder_sequence)

    def test_postorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Postorder Internal Node Iterator ###

    def test_postorder_internal_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_internal_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x]]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_internal_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_internal_node_iter(exclude_seed_node=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_internal_node_iter(exclude_seed_node=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Level-Order Node Iterator ###

    def test_levelorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.levelorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.levelorder_sequence)

    def test_levelorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.levelorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.levelorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### In-Order Node Iterator ###

    def test_inorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.inorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.inorder_sequence)

    def test_inorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.inorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.inorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Leaf Node Iterator ###

    def test_leaf_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.leaf_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.leaf_sequence)

    def test_leaf_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.leaf_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.leaf_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Age-Order Node Iterator ###

    def test_node_ages(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        tree.calc_node_ages()
        nodes = [nd for nd in tree.ageorder_node_iter()]
        for nd in nodes:
            self.assertEqual(nd.age, self.node_ages[nd.label])

    def test_ageorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.ageorder_sequence)

    def test_ageorder_node_iter_unfiltered_no_leaves(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter(include_leaves=False)]
        visited_labels = [nd.label for nd in nodes]
        expected = [label for label in self.ageorder_sequence if self.node_children[label]]
        self.assertSequenceEqual(visited_labels, expected)

    def test_ageorder_node_iter_unfiltered_reversed(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter(descending=True)]
        visited_labels = [nd.label for nd in nodes]
        nda = [ (self.node_ages[x], x) for x in self.preorder_sequence ]
        nda.sort(key=lambda x: x[0], reverse=True)
        exp = [x[1] for x in nda]
        self.assertSequenceEqual(visited_labels, exp)

    def test_leaf_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.ageorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.ageorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Preorder Edge Iterator ###

    def test_preorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.preorder_sequence)

    def test_preorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Preorder Internal Edge Iterator ###

    def test_preorder_internal_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x]]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(exclude_seed_edge=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(exclude_seed_edge=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Postorder Edge Iterator ###

    def test_postorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.postorder_sequence)

    def test_postorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Postorder Internal Edge Iterator ###

    def test_postorder_internal_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x]]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Level-Order Edge Iterator ###

    def test_levelorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.levelorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.levelorder_sequence)

    def test_levelorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.levelorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.levelorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### In-Order Edge Iterator ###

    def test_inorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.inorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.inorder_sequence)

    def test_inorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.inorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.inorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Leaf Edge Iterator ###

    def test_leaf_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.leaf_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertSequenceEqual(visited_labels, self.leaf_sequence)

    def test_leaf_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.leaf_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.leaf_sequence if self.node_edge_lengths[x] > 13]
        self.assertSequenceEqual(visited_labels, exp_labels)

    ### Special Iterators ###

    def test_child_node_iterator_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            expected_children = self.node_children[nd.label]
            children = [ch.label for ch in nd.child_node_iter()]
            self.assertSequenceEqual(children, expected_children)

    def test_child_node_iterator_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_children = [label for label in self.node_children[nd.label] if self.node_edge_lengths[label] > 13]
            children = [ch.label for ch in nd.child_node_iter(filter_fn=filter_fn)]
            self.assertEqual(children, expected_children)

    def test_child_edge_iterator_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            expected_children = self.node_children[nd.label]
            children = [edge.head_node.label for edge in nd.child_edge_iter()]
            self.assertSequenceEqual(children, expected_children)

    def test_child_edge_iterator_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.length > 13
        for nd in anodes:
            expected_children = [label for label in self.node_children[nd.label] if self.node_edge_lengths[label] > 13]
            children = [edge.head_node.label for edge in nd.child_edge_iter(filter_fn=filter_fn)]
            self.assertEqual(children, expected_children)

    def test_ancestor_iterator_exclusive_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=False)]
            expected_ancestors = self.node_ancestors[nd.label]
            self.assertSequenceEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_exclusive_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_ancestors = self.node_ancestors[nd.label]
            expected_ancestors = [nda for nda in expected_ancestors if self.node_edge_lengths[nda] > 13]
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=False, filter_fn=filter_fn)]
            self.assertEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_inclusive_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=True)]
            expected_ancestors = [nd.label] + list(self.node_ancestors[nd.label])
            self.assertEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_inclusive_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_ancestors = [nd.label] + list(self.node_ancestors[nd.label])
            expected_ancestors = [nda for nda in expected_ancestors if self.node_edge_lengths[nda] > 13]
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=True, filter_fn=filter_fn)]
            self.assertEqual(ancestors, expected_ancestors)

class TreeRootingState(dendropytest.ExtendedTestCase):

    def test_is_rooted(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_is_unrooted(self):
        self.assertFalse(self.fail_incomplete_tests())

class TestTreeApply(curated_test_tree.CuratedTestTree, unittest.TestCase):

    def test_apply(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        observed = []
        before_fn = lambda x: observed.append(("before", x.label))
        after_fn = lambda x: observed.append(("after", x.label))
        leaf_fn = lambda x: observed.append(("leaf", x.label))
        tree.apply(
                before_fn=before_fn,
                after_fn=after_fn,
                leaf_fn=leaf_fn)
        expected = [
            ("before", "a"),
            ("before", "b"),
            ("leaf", "i"),
            ("before", "e"),
            ("leaf", "j"),
            ("leaf", "k"),
            ("after", "e"),
            ("after", "b"),
            ("before", "c"),
            ("before", "g"),
            ("leaf", "l"),
            ("leaf", "m"),
            ("after", "g"),
            ("before", "f"),
            ("leaf", "n"),
            ("before", "h"),
            ("leaf", "o"),
            ("leaf", "p"),
            ("after", "h"),
            ("after", "f"),
            ("after", "c"),
            ("after", "a"),
                ]
        self.assertEqual(observed, expected)

if __name__ == "__main__":
    unittest.main()

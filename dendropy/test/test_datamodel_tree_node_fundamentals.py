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
Tests basic Node child management.
"""

import unittest
import dendropy
import copy
from dendropy.test.support import compare_and_validate

class TestNodeConstruction(unittest.TestCase):

    def test_basic_construction(self):
        taxon = dendropy.Taxon("z")
        nd = dendropy.Node(taxon=taxon, label="x", edge_length=1)
        self.assertIs(nd.taxon, taxon)
        self.assertEqual(nd.label, "x")
        edge = nd.edge
        self.assertEqual(edge.length, 1)
        self.assertIs(edge.head_node, nd)
        self.assertIs(edge.tail_node, None)

class NodeIdentity(unittest.TestCase):

    def setUp(self):
        taxon = dendropy.Taxon("a")
        self.n1 = dendropy.Node(label="a", taxon=taxon)
        self.n2 = dendropy.Node(label="a", taxon=taxon)

    def test_equal(self):
        # two distinct |Node| objects are never equal, even if all
        # member values are the same.
        self.assertNotEqual(self.n1, self.n2)
        self.assertIs(self.n1.taxon, self.n2.taxon)

    def test_hash_dict_membership(self):
        k = {}
        k[self.n1] = 1
        k[self.n2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.n1], 1)
        self.assertEqual(k[self.n2], 2)
        self.assertIn(self.n1, k)
        self.assertIn(self.n2, k)
        del k[self.n1]
        self.assertNotIn(self.n1, k)
        self.assertIn(self.n2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.n1: 1}
        k2 = {self.n2: 1}
        self.assertIn(self.n1, k1)
        self.assertIn(self.n2, k2)
        self.assertNotIn(self.n2, k1)
        self.assertNotIn(self.n1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.n1)
        k.add(self.n2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.n1, k)
        self.assertIn(self.n2, k)
        k.discard(self.n1)
        self.assertNotIn(self.n1, k)
        self.assertIn(self.n2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.n1: 1}
        k2 = {self.n2: 1}
        self.assertIn(self.n1, k1)
        self.assertIn(self.n2, k2)
        self.assertNotIn(self.n2, k1)
        self.assertNotIn(self.n1, k2)

class NodeCloning(compare_and_validate.Comparator, unittest.TestCase):

    def setUp(self):
        self.taxa = [dendropy.Taxon(label=label) for label in ["a", "b", "c", "d"]]
        self.n0 = dendropy.Node(label="0", taxon=self.taxa[0])
        self.c1 = dendropy.Node(label="1", taxon=None)
        self.c2 = dendropy.Node(label=None, taxon=self.taxa[1])
        self.c3 = dendropy.Node(label=None, taxon=None)
        self.c3 = dendropy.Node(label=None, taxon=self.taxa[2])
        self.p1 = dendropy.Node(label="-1", taxon=self.taxa[3])
        self.n0.parent_node = self.p1
        self.n0.set_child_nodes([self.c1, self.c2])
        self.c2.set_child_nodes([self.c3])
        self.nodes = [self.n0, self.c1, self.c2, self.c3, self.p1]
        for idx, nd in enumerate(self.nodes):
            if idx % 2 == 0:
                nd.edge.label = "E{}".format(idx)
                nd.edge.length = idx
            an1 = nd.annotations.add_new("a{}".format(idx),
                    "{}{}{}".format(nd.label, nd.taxon, idx))
            an2 = nd.annotations.add_bound_attribute("label")
            an3 = an1.annotations.add_bound_attribute("name")
            ae1 = nd.edge.annotations.add_new("a{}".format(idx),
                    "{}{}".format(nd.edge.label, idx))
            ae2 = nd.edge.annotations.add_bound_attribute("label")
            ae3 = ae1.annotations.add_bound_attribute("name")

    def test_unsupported_copy(self):
        with self.assertRaises(TypeError):
            self.n0.clone(0)
        with self.assertRaises(TypeError):
            copy.copy(self.n0)
        with self.assertRaises(TypeError):
            self.n0.clone(1)
        with self.assertRaises(TypeError):
            self.n0.taxon_namespace_scoped_copy()

    def test_deepcopy(self):
        for clone in (
                self.n0.clone(2),
                copy.deepcopy(self.n0),
                ):
            self.compare_distinct_nodes(
                    clone, self.n0,
                    taxon_namespace_scoped=False,
                    compare_tree_annotations=True)

class TestNodeSetChildNodes(unittest.TestCase):

    def test_set_child_nodes(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        for nd in assigned_ch:
            x = [dendropy.Node(label=c) for c in ["s1", "s2"]]
            nd.set_child_nodes(x)
            nd._expected_children = x
        parent.set_child_nodes(assigned_ch)
        for ch in parent._child_nodes:
            self.assertIn(ch, assigned_ch)
            self.assertIs(ch._parent_node, parent)
            self.assertIs(ch.edge.tail_node, parent)
            self.assertIs(ch.edge.head_node, ch)
            self.assertEqual(len(ch._child_nodes), len(ch._expected_children))
            for sch in ch._child_nodes:
                self.assertIn(sch, ch._expected_children)
                self.assertIs(sch._parent_node, ch)
                self.assertIs(sch.edge.tail_node, ch)
                self.assertIs(sch.edge.head_node, sch)
        for ch in assigned_ch:
            self.assertTrue(ch in parent._child_nodes)

    def test_add_child(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        for nd in assigned_ch:
            x = [dendropy.Node(label=c) for c in ["s1", "s2"]]
            for y in x:
                nd.add_child(y)
            nd._expected_children = x
        for ch in assigned_ch:
            k = parent.add_child(ch)
            self.assertIs(k, ch)
        for ch in parent._child_nodes:
            self.assertIn(ch, assigned_ch)
            self.assertIs(ch._parent_node, parent)
            self.assertIs(ch.edge.tail_node, parent)
            self.assertIs(ch.edge.head_node, ch)
            self.assertEqual(len(ch._child_nodes), len(ch._expected_children))
            for sch in ch._child_nodes:
                self.assertIn(sch, ch._expected_children)
                self.assertIs(sch._parent_node, ch)
                self.assertIs(sch.edge.tail_node, ch)
                self.assertIs(sch.edge.head_node, sch)
        for ch in assigned_ch:
            self.assertTrue(ch in parent._child_nodes)

    def test_insert_child_at_pos(self):
        new_child_labels = ["c1", "c2", "c3"]
        insert_ch_label = "x1"
        for pos in range(len(new_child_labels)+1):
            parent = dendropy.Node(label="parent")
            assigned_ch = [dendropy.Node(label=c) for c in new_child_labels]
            parent.set_child_nodes(assigned_ch)
            insert_ch = dendropy.Node(label=insert_ch_label)
            parent.insert_child(pos, insert_ch)
            x = 0
            for idx, ch in enumerate(parent._child_nodes):
                if idx == pos:
                    self.assertEqual(ch.label, insert_ch_label)
                else:
                    self.assertEqual(ch.label, new_child_labels[x])
                    x += 1

    def test_redundant_set(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        parent.set_child_nodes(assigned_ch)
        parent.set_child_nodes(assigned_ch)
        ch2 = assigned_ch + assigned_ch
        parent.set_child_nodes(ch2)
        self.assertEqual(parent._child_nodes, assigned_ch)
        for nd in parent.child_node_iter():
            self.assertIs(nd.parent_node, parent)

    def test_redundant_insert_child_at_pos(self):
        new_child_labels = ["c1", "c2", "c3"]
        for child_to_insert_idx in range(len(new_child_labels)):
            for insertion_idx in range(len(new_child_labels)):
                parent = dendropy.Node(label="parent")
                assigned_ch = [dendropy.Node(label=c) for c in new_child_labels]
                parent.set_child_nodes(assigned_ch)
                self.assertEqual(parent._child_nodes, assigned_ch)
                insert_ch = assigned_ch[child_to_insert_idx]
                parent.insert_child(insertion_idx, insert_ch)
                self.assertEqual(len(parent._child_nodes), len(assigned_ch))
                self.assertEqual(len(set(parent._child_nodes)), len(parent._child_nodes))
                x = 0
                for idx, ch in enumerate(parent._child_nodes):
                    if idx == insertion_idx:
                        self.assertIs(ch, insert_ch)
                    self.assertIn(ch, assigned_ch)
                for ch in assigned_ch:
                    self.assertIn(ch, parent._child_nodes)
                    self.assertEqual(parent._child_nodes.count(ch), 1)

    def test_redundant_add(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        parent.set_child_nodes(assigned_ch)
        for ch in assigned_ch:
            parent.add_child(ch)
            self.assertEqual(parent._child_nodes, assigned_ch)
            for nd in parent.child_node_iter():
                self.assertIs(nd.parent_node, parent)

    def test_parent_node_setting(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        for ch in assigned_ch:
            ch.parent_node = parent
        for ch in assigned_ch:
            self.assertEqual(parent._child_nodes, assigned_ch)
            for nd in parent.child_node_iter():
                self.assertIs(nd.parent_node, parent)
                self.assertIs(nd.edge.tail_node, parent)
                self.assertIs(nd.edge.head_node, nd)

    def test_new_child(self):
        parent = dendropy.Node(label="parent")
        new_child_labels = ["c1", "c2", "c3"]
        sub_child_labels = ["s1", "s2"]
        for label in new_child_labels:
            nd = parent.new_child(label=label)
            for y in sub_child_labels:
                x = nd.new_child(label=y)
                self.assertTrue(isinstance(x, dendropy.Node))
                self.assertEqual(x.label, y)
                self.assertIs(x.parent_node, nd)
                self.assertIs(x.edge.head_node, x)
                self.assertIs(x.edge.tail_node, nd)
        self.assertEqual(len(parent._child_nodes), len(new_child_labels))
        for ch in parent._child_nodes:
            self.assertIn(ch.label, new_child_labels)
            self.assertIs(ch._parent_node, parent)
            self.assertIs(ch.edge.tail_node, parent)
            self.assertIs(ch.edge.head_node, ch)
            for sch in ch._child_nodes:
                self.assertIn(sch.label, sub_child_labels)
                self.assertIs(sch._parent_node, ch)
                self.assertIs(sch.edge.tail_node, ch)
                self.assertIs(sch.edge.head_node, sch)

    def test_new_child_at_pos(self):
        new_child_labels = ["c1", "c2", "c3"]
        insert_ch_label = "x1"
        for pos in range(len(new_child_labels)+1):
            parent = dendropy.Node(label="parent")
            assigned_ch = [dendropy.Node(label=c) for c in new_child_labels]
            parent.set_child_nodes(assigned_ch)
            parent.insert_new_child(pos, label=insert_ch_label)
            x = 0
            for idx, ch in enumerate(parent._child_nodes):
                if idx == pos:
                    self.assertEqual(ch.label, insert_ch_label)
                else:
                    self.assertEqual(ch.label, new_child_labels[x])
                    x += 1

    def test_remove_child(self):
        assigned_child_labels = ["c1", "c2", "c3"]
        for remove_idx in range(len(assigned_child_labels)):
            parent = dendropy.Node(label="parent")
            assigned_ch = [dendropy.Node(label=c) for c in assigned_child_labels]
            parent.set_child_nodes(assigned_ch)
            ch_nodes = list(parent._child_nodes)
            to_remove = ch_nodes[remove_idx]
            x = parent.remove_child(to_remove)
            self.assertIs(to_remove, x)
            ch_nodes.remove(to_remove)
            ch_nodes2 = list(parent._child_nodes)
            self.assertEqual(ch_nodes, ch_nodes2)

    def test_edge_head_node_setting(self):
        node = dendropy.Node(label="x")
        edge1 = node.edge
        edge2 = dendropy.Edge()
        node.edge = edge2
        self.assertIs(node.edge, edge2)
        self.assertIs(node.edge.head_node, node)

if __name__ == "__main__":
    unittest.main()

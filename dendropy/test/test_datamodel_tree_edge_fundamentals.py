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
Tests basic Edge operations.
"""

import copy
import unittest
from dendropy.test.support import compare_and_validate
import dendropy

class EdgeIdentity(unittest.TestCase):

    def setUp(self):
        self.e1 = dendropy.Edge(label="a")
        self.e2 = dendropy.Edge(label="a")

    def test_equal(self):
        # two distinct |Edge| objects are never equal, even if all
        # member values are the same.
        self.assertNotEqual(self.e1, self.e2)

    def test_hash_dict_membership(self):
        k = {}
        k[self.e1] = 1
        k[self.e2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.e1], 1)
        self.assertEqual(k[self.e2], 2)
        self.assertIn(self.e1, k)
        self.assertIn(self.e2, k)
        del k[self.e1]
        self.assertNotIn(self.e1, k)
        self.assertIn(self.e2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.e1: 1}
        k2 = {self.e2: 1}
        self.assertIn(self.e1, k1)
        self.assertIn(self.e2, k2)
        self.assertNotIn(self.e2, k1)
        self.assertNotIn(self.e1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.e1)
        k.add(self.e2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.e1, k)
        self.assertIn(self.e2, k)
        k.discard(self.e1)
        self.assertNotIn(self.e1, k)
        self.assertIn(self.e2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.e1: 1}
        k2 = {self.e2: 1}
        self.assertIn(self.e1, k1)
        self.assertIn(self.e2, k2)
        self.assertNotIn(self.e2, k1)
        self.assertNotIn(self.e1, k2)

class EdgeCloning(compare_and_validate.Comparator, unittest.TestCase):

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
        self.e0 = self.n0._edge

    def test_unsupported_copy(self):
        with self.assertRaises(TypeError):
            self.e0.clone(0)
        with self.assertRaises(TypeError):
            copy.copy(self.e0)
        with self.assertRaises(TypeError):
            self.e0.clone(1)
        with self.assertRaises(TypeError):
            self.e0.taxon_namespace_scoped_copy()

    def test_deepcopy(self):
        for clone in (
                self.e0.clone(2),
                copy.deepcopy(self.e0),
                ):
            self.compare_distinct_nodes(
                    clone._head_node, self.n0,
                    taxon_namespace_scoped=False,
                    compare_tree_annotations=True)

class EdgeNodeManagement(unittest.TestCase):


    def test_edge_tail_node_setting(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(label=c) for c in ["c1", "c2", "c3"]]
        for ch in assigned_ch:
            ch.edge.tail_node = parent
        for ch in assigned_ch:
            self.assertEqual(parent._child_nodes, assigned_ch)
            for nd in parent.child_node_iter():
                self.assertIs(nd.parent_node, parent)
                self.assertIs(nd.edge.tail_node, parent)
                self.assertIs(nd.edge.head_node, nd)

    def test_edge_head_node_setting(self):

        node1 = dendropy.Node()
        parent_node1 = dendropy.Node()
        node1.parent_node = parent_node1
        self.assertIs(node1._parent_node, parent_node1)
        edge1 = node1.edge
        self.assertIs(edge1.head_node, node1)
        self.assertIs(edge1.tail_node, parent_node1)

        node2 = dendropy.Node()
        parent_node2 = dendropy.Node()
        node2.parent_node = parent_node2
        edge2 = node2.edge
        self.assertIs(edge2.head_node, node2)
        self.assertIs(edge2.tail_node, parent_node2)

        new_edge1 = dendropy.Edge()
        new_edge1.head_node = node1
        self.assertIs(node1.edge, new_edge1)
        self.assertIs(new_edge1.head_node, node1)
        self.assertIs(new_edge1.tail_node, parent_node1)

        edge2.head_node = node1
        self.assertIs(node1.edge, edge2)
        self.assertIs(edge2.head_node, node1)
        self.assertIs(edge2.tail_node, parent_node1)

if __name__ == "__main__":
    unittest.main()

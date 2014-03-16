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
Tests basic Node operations.
"""

import unittest
import dendropy

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

class TestNodeSetChildNodes(unittest.TestCase):

    def test_set_child_node(self):
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
            parent.add_child(ch)
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

    def test_add_child_at_pos(self):
        new_child_labels = ["c1", "c2", "c3"]
        insert_ch_label = "x1"
        for pos in range(len(new_child_labels)+1):
            parent = dendropy.Node(label="parent")
            assigned_ch = [dendropy.Node(label=c) for c in new_child_labels]
            parent.set_child_nodes(assigned_ch)
            insert_ch = dendropy.Node(label=insert_ch_label)
            parent.add_child(insert_ch, pos)
            x = 0
            for idx, ch in enumerate(parent._child_nodes):
                if idx == pos:
                    self.assertEqual(ch.label, insert_ch_label)
                else:
                    self.assertEqual(ch.label, new_child_labels[x])
                    x += 1

    def test_new_child(self):
        parent = dendropy.Node(label="parent")
        new_child_labels = ["c1", "c2", "c3"]
        sub_child_labels = ["s1", "s2"]
        for label in new_child_labels:
            nd = parent.new_child(label=label)
            for y in sub_child_labels:
                nd.new_child(label=y)
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

if __name__ == "__main__":
    unittest.main()

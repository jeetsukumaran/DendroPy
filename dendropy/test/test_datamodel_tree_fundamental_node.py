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
        assigned_ch = [dendropy.Node(c) for c in ["c1", "c2", "c3"]]
        for nd in assigned_ch:
            x = [dendropy.Node(c) for c in ["s1", "s2"]]
            nd.set_child_nodes(x)
            nd._expected_children = x
        parent.set_child_nodes(assigned_ch)
        for ch in parent._child_nodes:
            self.assertIn(ch, assigned_ch)
            self.assertIs(ch._parent_node, parent)
            self.assertIs(ch.edge.tail_node, parent)
            self.assertIs(ch.edge.head_node, ch)
            for sch in ch._child_nodes:
                self.assertIn(sch, ch._expected_children)
                self.assertIs(sch._parent_node, ch)
                self.assertIs(sch.edge.tail_node, ch)
                self.assertIs(sch.edge.head_node, sch)
        for ch in assigned_ch:
            self.assertTrue(ch in parent._child_nodes)

    def test_add_child(self):
        parent = dendropy.Node(label="parent")
        assigned_ch = [dendropy.Node(c) for c in ["c1", "c2", "c3"]]
        for nd in assigned_ch:
            x = [dendropy.Node(c) for c in ["s1", "s2"]]
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
            for sch in ch._child_nodes:
                self.assertIn(sch, ch._expected_children)
                self.assertIs(sch._parent_node, ch)
                self.assertIs(sch.edge.tail_node, ch)
                self.assertIs(sch.edge.head_node, sch)
        for ch in assigned_ch:
            self.assertTrue(ch in parent._child_nodes)

if __name__ == "__main__":
    unittest.main()

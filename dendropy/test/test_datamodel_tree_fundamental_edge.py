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
Tests basic Edge operations.
"""

import copy
import unittest
from dendropy.test.support import compare_and_validate
import dendropy

class EdgeIdentity(unittest.TestCase):

    def setUp(self):
        self.e1 = dendropy.Edge("a")
        self.e2 = dendropy.Edge("a")

    def test_equal(self):
        # two distinct :class:`Edge` objects are never equal, even if all
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

    def test_hash_set_membership(self):
        k = set()
        k.add(self.e1)
        k.add(self.e2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.e1, k)
        self.assertIn(self.e2, k)

# note that compare_and_validate.AnnotableComparator must be listed first,
# otherwise setUp will not be called
class EdgeCloning(compare_and_validate.AnnotableComparator, unittest.TestCase):

    def setUp(self):
        self.e0 = dendropy.Edge(label="0", length=1)
        self.e0.tail_node = dendropy.Node("-1")
        self.e0.head_node = dendropy.Node("0")
        self.e0.rootedge = True

    def test_copy_from_another(self):
        edge_copy = self.e0.clone(0)
        self.assertIsNot(self.e0, edge_copy)
        self.assertNotEqual(self.e0, edge_copy)
        self.assertEqual(self.e0.label, edge_copy.label)
        self.assertEqual(self.e0.length, edge_copy.length)
        self.assertEqual(self.e0.rootedge, edge_copy.rootedge)
        self.assertIs(edge_copy.head_node, self.e0.head_node)
        self.assertIs(edge_copy.tail_node, self.e0.tail_node)

    def test_copy_from_another_with_simple_annotations(self):
        self.e0 = dendropy.Edge("a")
        self.e0.annotations.add_new("a", 0)
        self.e0.annotations.add_new("b", 1)
        self.e0.annotations.add_new("c", 3)
        edge_copy = self.e0.clone(0)
        self.assertIsNot(self.e0, edge_copy)
        self.assertNotEqual(self.e0, edge_copy)
        self.assertEqual(self.e0.label, edge_copy.label)
        self.assertTrue(hasattr(self.e0, "annotations"))
        self.assertTrue(hasattr(edge_copy, "annotations"))
        self.assertEqual(len(self.e0.annotations), len(edge_copy.annotations))
        self.compare_annotables(self.e0, edge_copy)

    def test_copy_from_another_with_complex_annotations(self):
        self.e0 = dendropy.Edge("a")
        self.e0.annotations.add_new("a", 0)
        b = self.e0.annotations.add_new("b", (self.e0, "label"), is_attribute=True)
        b.annotations.add_new("c", 3)
        edge_copy = self.e0.clone(0)
        self.assertIsNot(self.e0, edge_copy)
        self.assertNotEqual(self.e0, edge_copy)
        self.assertEqual(self.e0.label, edge_copy.label)
        self.assertTrue(hasattr(self.e0, "annotations"))
        self.assertTrue(hasattr(edge_copy, "annotations"))
        self.assertEqual(len(self.e0.annotations), len(edge_copy.annotations))
        self.compare_annotables(self.e0, edge_copy)
        self.e0.label = "x"
        edge_copy.label = "y"
        self.assertEqual(self.e0.annotations[1].value, "x")
        self.assertEqual(edge_copy.annotations[1].value, "y")

if __name__ == "__main__":
    unittest.main()

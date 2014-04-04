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

    def test_deepcopy_from_another(self):
        e1 = dendropy.Edge("a")
        e2 = copy.deepcopy(e1)
        self.assertIsNot(e1, e2)
        self.assertNotEqual(e1, e2)
        self.assertEqual(e1.label, e2.label)

    def test_deepcopy_from_another_with_simple_annotations(self):
        e1 = dendropy.Edge("a")
        e1.annotations.add_new("a", 0)
        e1.annotations.add_new("b", 1)
        e1.annotations.add_new("c", 3)
        e2 = copy.deepcopy(e1)
        self.assertIsNot(e1, e2)
        self.assertNotEqual(e1, e2)
        self.assertEqual(e1.label, e2.label)
        self.assertTrue(hasattr(e1, "annotations"))
        self.assertTrue(hasattr(e2, "annotations"))
        self.assertEqual(len(e1.annotations), len(e2.annotations))
        self.compare_annotables(e1, e2)

    def test_deepcopy_from_another_with_complex_annotations(self):
        e1 = dendropy.Edge("a")
        e1.annotations.add_new("a", 0)
        b = e1.annotations.add_new("b", (e1, "label"), is_attribute=True)
        b.annotations.add_new("c", 3)
        e2 = copy.deepcopy(e1)
        self.assertIsNot(e1, e2)
        self.assertNotEqual(e1, e2)
        self.assertEqual(e1.label, e2.label)
        self.assertTrue(hasattr(e1, "annotations"))
        self.assertTrue(hasattr(e2, "annotations"))
        self.assertEqual(len(e1.annotations), len(e2.annotations))
        self.compare_annotables(e1, e2)
        e1.label = "x"
        e2.label = "y"
        self.assertEqual(e1.annotations[1].value, "x")
        self.assertEqual(e2.annotations[1].value, "y")

if __name__ == "__main__":
    unittest.main()

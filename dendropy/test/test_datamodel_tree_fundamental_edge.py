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
        self.e1 = dendropy.Edge(label="a")
        self.e2 = dendropy.Edge(label="a")

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

if __name__ == "__main__":
    unittest.main()

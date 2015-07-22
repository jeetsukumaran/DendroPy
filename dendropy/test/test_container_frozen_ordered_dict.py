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
ReadonlyOrdedr tests.
"""

import copy
import collections
import unittest
from dendropy.utility import container

class FrozenOrderedDictTest(unittest.TestCase):

    def setUp(self):
        self.src = collections.OrderedDict()
        self.src["a"] = 1
        self.src["b"] = 2
        self.src["c"] = 3
        self.d = container.FrozenOrderedDict(self.src)

    def test_basic_readonly_iteration_and_access(self):
        for key1, key2 in zip(self.src, self.d):
            self.assertEqual(key1, key2)
            key = key1
            self.assertTrue(key in self.d)
            self.assertEqual(self.d[key], self.src[key])

    def test_keys(self):
        k1 = list(self.src.keys())
        k2 = list(self.d.keys())
        self.assertEqual(k1, k2)

    def test_values(self):
        v1 = list(self.src.values())
        v2 = list(self.d.values())
        self.assertEqual(v1, v2)

    def test_readonly(self):
        k = list(self.d.keys())[0]
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            self.d[k] = 1
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            del self.d[k]
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            self.d.pop(k)
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            self.d.clear()
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            self.d.update({})
        with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
            self.d.fromkeys([1,2,3])

    def test_deepcopy(self):
        d2 = copy.deepcopy(self.d)
        self.assertIsNot(d2, self.d)
        for key1, key2 in zip(d2, self.d):
            self.assertEqual(key1, key2)
            key = key1
            self.assertTrue(key in self.d)
            self.assertEqual(self.d[key], d2[key])

    def test_copy(self):
        d2 = copy.copy(self.d)
        self.assertEqual(d2, self.d)

    def test_copy_construction(self):
        d2 = collections.OrderedDict(self.d)
        self.assertEqual(d2, self.d)

if __name__ == "__main__":
    unittest.main()

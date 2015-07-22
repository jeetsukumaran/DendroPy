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
OrderedSet tests.
"""

import random
import copy
import unittest
from dendropy.utility import container

class A(object):
    def __hash__(self):
        return id(self)
    def __eq__(self, other):
        return self is other

class B(object):
    def __init__(self, v):
        self.data = [v * random.randint(-1000, 100) for i in range(10)]
    def __hash__(self):
        return id(self)
    def __eq__(self, other):
        return self.data == other.data

class TestOrderedSetIdentity(unittest.TestCase):

    def setUp(self):
        self.t1 = container.OrderedSet(["a", "b"])
        self.t2 = container.OrderedSet(["a", "b"])

    def test_hash_dict_membership(self):
        k = {}
        k[self.t1] = 1
        k[self.t2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.t1], 1)
        self.assertEqual(k[self.t2], 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        del k[self.t1]
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.t1)
        k.add(self.t2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        k.discard(self.t1)
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

class TestOrderedSetCollectionsManagement(unittest.TestCase):

    def test_basic_adding(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        ordered_set = container.OrderedSet()
        for item in items:
            ordered_set.add(item)
        for x in range(100):
            ordered_set.add( random.choice(items) )
        result = [item for item in ordered_set]
        self.assertEqual(result, items)
        self.assertEqual(len(result), len(items))

    def test_constructor(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        ordered_set = container.OrderedSet(items)
        for x in range(100):
            ordered_set.add( random.choice(items) )
        result = [item for item in ordered_set]
        self.assertEqual(result, items)
        self.assertEqual(len(result), len(items))

    def test_basic_getting(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        ordered_set = container.OrderedSet()
        for item in items:
            ordered_set.add(item)
        for idx in range(len(items)):
            self.assertIs(ordered_set[idx], items[idx])

    def test_del(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        for x in range(len(items)):
            ordered_set = container.OrderedSet(items)
            zz = list(items)
            while len(zz):
                i = random.randint(0, len(zz)-1)
                self.assertIs(ordered_set[i], zz[i])
                del ordered_set[i]
                del zz[i]
                result = [item for item in ordered_set]
                self.assertEqual(result, zz)
                self.assertEqual(len(result), len(zz))

    def test_discard(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        for x in range(len(items)):
            ordered_set = container.OrderedSet(items)
            zz = list(items)
            while len(zz):
                k = random.choice(zz)
                ordered_set.discard(k)
                zz.remove(k)
                result = [item for item in ordered_set]
                self.assertEqual(result, zz)
                self.assertEqual(len(result), len(zz))

    def test_pop_back(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        ordered_set = container.OrderedSet(items)
        while ordered_set:
            item1 = ordered_set.pop()
            item2 = items.pop(-1)
            self.assertIs(item1, item2)
            self.assertEqual(len(ordered_set), len(items))
            result = [item for item in ordered_set]
            self.assertEqual(result, items)
        self.assertEqual(len(ordered_set), 0)

    def test_pop_front(self):
        items = [42, 3.14, "hello", object(), A(), frozenset([1,2,3]), A()]
        ordered_set = container.OrderedSet(items)
        while ordered_set:
            item1 = ordered_set.pop(False)
            item2 = items.pop(0)
            self.assertIs(item1, item2)
            self.assertEqual(len(ordered_set), len(items))
            result = [item for item in ordered_set]
            self.assertEqual(result, items)
        self.assertEqual(len(ordered_set), 0)

class TestOrderedSetDeepCopy(unittest.TestCase):

    def setUp(self):
        self.s1 = container.OrderedSet([0,1,2,3,"a","b","c",])
        items = [B(1), B(0), B(2), frozenset([1,2,3]), B(3),]
        self.c1 = container.OrderedSet(items)

    def test_simple_element_deepcopy(self):
        s2 = copy.deepcopy(self.s1)
        self.assertEqual(s2, self.s1)

    def test_deepcopy(self):
        c2 = copy.deepcopy(self.c1)
        self.assertEqual(len(c2), len(self.c1))
        for item2, item1 in zip(c2, self.c1):
            self.assertIsNot(item2, item1)
            self.assertEqual(item2, item1)
        for k in self.c1.__dict__:
            self.assertIn(k, c2.__dict__)
            self.assertIsNot(c2.__dict__[k], self.c1.__dict__[k])

if __name__ == "__main__":
    unittest.main()

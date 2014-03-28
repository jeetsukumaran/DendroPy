#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
OrderedSet tests.
"""

import random
import unittest
from dendropy.utility import container

class A(object):
    pass

class TestOrderedSet(unittest.TestCase):

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

if __name__ == "__main__":
    unittest.main()

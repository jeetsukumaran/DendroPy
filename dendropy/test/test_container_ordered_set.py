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

if __name__ == "__main__":
    unittest.main()

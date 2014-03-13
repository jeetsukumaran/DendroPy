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
Annotation tests.
"""

import collections
import unittest
import copy
from dendropy.datamodel import base

class TestObject(base.Annotable):
    pass

class DummyX(TestObject):
    def __init__(self, data=None):
        self.data = data

class AnnotableDeepCopyTester(unittest.TestCase):

    def compare(self, x1, x2):
        self.assertIsNot(x1, x2)
        if not hasattr(x1, "_annotations"):
            self.assertTrue(not hasattr(x2, "_annotations"))
            return
        self.assertIs(x1._annotations.target, x1)
        self.assertIs(x2._annotations.target, x2)
        self.assertTrue(hasattr(x2, "_annotations"))
        self.assertIsNot(x1._annotations, x2._annotations)
        self.assertEqual(len(x1._annotations), len(x2._annotations))
        for a1, a2 in zip(x1._annotations, x2._annotations):
            self.assertIsNot(a1, a2)
            for k in a1.__dict__:
                self.assertIn(k, a2.__dict__)
                v1 = a1.__dict__[k]
                v2 = a2.__dict__[k]
                if isinstance(v1, TestObject):
                    self.assertNotIs(v1, v2)
                else:
                    self.assertEqual(v1, v2)
                self.compare(a1, a2)

    def test_deep_copy(self):
        x1 = DummyX()
        a1 = x1.annotations.add_new(name="a1", value="1")
        a2 = x1.annotations.add_new(name="a2", value="2")
        a3 = x1.annotations.add_new(name="a3", value="3")

        x2 = copy.deepcopy(x1)

        self.compare(x1, x2)


if __name__ == "__main__":
    unittest.main()

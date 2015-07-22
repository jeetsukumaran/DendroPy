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
Annotation tests.
"""

import collections
import unittest
import copy
from dendropy.datamodel import basemodel
from dendropy.test.support import compare_and_validate

class TestObject(basemodel.Annotable, basemodel.DataObject):
    pass

class DummyX(TestObject):
    def __init__(self, data=None):
        self.data = data

class AnnotableDeepCopyTester(compare_and_validate.Comparator, unittest.TestCase):

    def test_deep_copy(self):
        x1 = DummyX()
        a1 = x1.annotations.add_new(name="a1", value="1")
        a2 = x1.annotations.add_new(name="a2", value="2")
        a3 = x1.annotations.add_new(name="a3", value="3")
        x2 = copy.deepcopy(x1)
        self.compare_distinct_annotables(x1, x2)

    def test_nested_deep_copy(self):
        x1 = DummyX()
        a1 = x1.annotations.add_new(name="a1", value="1")
        a2 = x1.annotations.add_new(name="a2", value="2")
        a3 = x1.annotations.add_new(name="a3", value="3")
        a4 = a3.annotations.add_new(name="a4", value="4")
        a5 = a3.annotations.add_new(name="a5", value="5")
        a6 = a5.annotations.add_new(name="a6", value="6")
        x2 = copy.deepcopy(x1)
        self.compare_distinct_annotables(x1, x2)

if __name__ == "__main__":
    unittest.main()

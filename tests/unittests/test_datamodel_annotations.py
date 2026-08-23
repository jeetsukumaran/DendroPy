#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Annotation tests.
"""

import collections
import unittest
import copy
from dendropy.datamodel import basemodel
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import compare_and_validate

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

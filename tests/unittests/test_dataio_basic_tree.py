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

import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import mockreader
from dendropy import dataio
import dendropy

import unittest

class DendropyDataIOTestMockTreeReader(mockreader.MockReader):

    def __init__(self):
        mockreader.MockReader.__init__(self)

    def process_read_call(self):
        tns = self.taxon_namespace_factory(label="test1")
        tree_list = self.tree_list_factory(label="test1", taxon_namespace=tns)
        tree = tree_list.new_tree()
        product = self.Product(
                taxon_namespaces=[tns],
                tree_lists=[tree_list],
                char_matrices=None)
        return product

class MockTestTreeTypeDerivedFromDendropyTree(dendropy.Tree):
    pass


class TestCustomTreeType(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        dataio.register_reader("dendropy_test_mock_tree_reader", DendropyDataIOTestMockTreeReader)

    def test_get_from(self):
        tree = MockTestTreeTypeDerivedFromDendropyTree.get_from_string("", "dendropy_test_mock_tree_reader")
        self.assertEqual(type(tree), MockTestTreeTypeDerivedFromDendropyTree)

if __name__ == "__main__":
    unittest.main()


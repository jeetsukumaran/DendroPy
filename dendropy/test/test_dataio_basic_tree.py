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

from dendropy.test.support import mockreader
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



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


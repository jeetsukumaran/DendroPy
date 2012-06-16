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
Test cloning of dataobjects.
"""

import unittest
import copy

from dendropy.test.support import pathmap
from dendropy.test.support import datatest
import dendropy


class TestTreeCloning(datatest.AnnotatedDataObjectVerificationTestCase):

    def setUp(self):
        s = pathmap.tree_source_stream("pythonidae.annotated.nexml")
        self.dataset = dendropy.DataSet.get_from_stream(s, "nexml")

    def testDeepCopy(self):
        tree1 = self.dataset.tree_lists[0][0]
        tree2 = copy.deepcopy(tree1)
        self.assertDistinctButEqualTree(tree1, tree2, distinct_taxa=False)

    def testCopyConstruction(self):
        tree1 = self.dataset.tree_lists[0][0]
        tree2 = dendropy.Tree(tree1)
        self.assertDistinctButEqualTree(tree1, tree2, distinct_taxa=False)

class TestTreeListCloning(datatest.AnnotatedDataObjectVerificationTestCase):

    def setUp(self):
        s = pathmap.tree_source_stream("pythonidae.annotated.nexml")
        self.dataset = dendropy.DataSet.get_from_stream(s, "nexml")

    def testDeepCopy(self):
        tree_list1 = self.dataset.tree_lists[0]
        tree_list2 = copy.deepcopy(tree_list1)
        self.assertDistinctButEqualTreeList(tree_list1, tree_list2, distinct_taxa=False)

    def testCopyConstruction(self):
        tree_list1 = self.dataset.tree_lists[0]
        tree_list2 = dendropy.TreeList(tree_list1)
        self.assertDistinctButEqualTreeList(tree_list1, tree_list2, distinct_taxa=False)


if __name__ == "__main__":
    unittest.main()

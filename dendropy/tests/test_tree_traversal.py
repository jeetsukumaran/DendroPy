#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Tests composition and traversal of trees.
"""

import unittest
from dendropy.utility import messaging
import dendropy

_LOG = messaging.get_logger(__name__)

class TreeTraversalChecker(unittest.TestCase):

    def setUp(self):
        """
        Sets up tree:
            ((((T1, ((T2, T3 )i5, T4 )i4 )i3, (T5, T6, T13 )i6 )i2, ((T7, (T8, T9 )i9 )i8, T10 )i7 )i1, T14, (T11, T12 )i10 )i0;
        """
        self.tree = dendropy.Tree(oid='tree1')
        self.tree.seed_node.oid = 'i0'
        node_i1 = self.tree.seed_node.new_child(oid='i1')
        node_i2 = node_i1.new_child(oid='i2')
        node_i3 = node_i2.new_child(oid='i3')
        node_i3.new_child(oid='T1')
        node_i4 = node_i3.new_child(oid='i4')
        node_i5 = node_i4.new_child(oid='i5')
        node_i5.new_child(oid='T2')
        node_i5.new_child(oid='T3')
        node_i4.new_child(oid='T4')
        node_i6 = node_i2.new_child(oid='i6')
        node_i6.new_child(oid='T5')
        node_i6.new_child(oid='T6')
        node_i7 = node_i1.new_child(oid='i7')
        node_i8 = node_i7.new_child(oid='i8')
        node_i8.new_child(oid='T7')
        node_i9 = node_i8.new_child(oid='i9')
        node_i9.new_child(oid='T8')
        node_i9.new_child(oid='T9')
        node_i7.new_child(oid='T10')
        self.tree.seed_node.new_child(oid='T14')
        node_i10 = self.tree.seed_node.new_child(oid='i10')
        node_i10.new_child(oid='T11')
        node_i10.new_child(oid='T12')
        node_i6.new_child(oid='T13')
        self.tree.debug_check_tree(_LOG)

class TestLeafIter(TreeTraversalChecker):

    def runTest(self):
        result = [node.oid for node in self.tree.leaf_iter()]
        expected = ["T1", "T2", "T3", "T4", "T5", "T6", "T13", "T7", "T8", "T9", "T10", "T14", "T11", "T12"]
        self.assertEqual(result, expected)

class TestPreorderIter(TreeTraversalChecker):

    def runTest(self):
        result = [node.oid for node in self.tree.preorder_node_iter()]
        expected = ["i0", "i1", "i2", "i3", "T1", "i4", "i5", "T2", "T3", "T4", "i6", "T5", "T6", "T13", "i7", "i8", "T7", "i9", "T8", "T9", "T10", "T14", "i10", "T11", "T12"]
        self.assertEqual(result, expected)

class TestPostorderIter(TreeTraversalChecker):

    def runTest(self):
        result = [node.oid for node in self.tree.postorder_node_iter()]
        expected = ["T1", "T2", "T3", "i5", "T4", "i4", "i3", "T5", "T6", "T13", "i6", "i2", "T7", "T8", "T9", "i9", "i8", "T10", "i7", "i1", "T14", "T11", "T12", "i10", "i0"]
        self.assertEqual(result, expected)

class TestLevelorderIter(TreeTraversalChecker):

    def runTest(self):
        result = [node.oid for node in self.tree.level_order_node_iter()]
        expected = ["i0", "i1", "T14", "i10", "i2", "i7", "T11", "T12", "i3", "i6", "i8", "T10", "T1", "i4", "T5", "T6", "T13", "T7", "i9", "i5", "T4", "T8", "T9", "T2", "T3"]
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()

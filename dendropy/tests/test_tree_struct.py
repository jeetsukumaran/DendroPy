#! /usr/bin/env python

############################################################################
##  test_tree_struct.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
Tests tree traversal.
"""

import unittest
from dendropy import get_logger
import dendropy.tests
_LOG = get_logger("Tree Structure and Traversal")

### MODULE THAT WE ARE TESTING ###
from dendropy import trees
### MODULE THAT WE ARE TESTING ###

class TreeIterTest(unittest.TestCase):

    def build_fixed_basic_tree1(self):
        """
        Returns a bifurcating basic tree looking like:

        http://en.wikipedia.org/wiki/Tree_traversal
        (   (   A,
                (   C,
                    E
                )D
            )B,
            (   (   H
                )I
            )G
        )F
        Preorder traversal sequence:    F, B, A, D, C, E, G, I, H
        Inorder traversal sequence:     A, B, C, D, E, F, G, H, I
        Postorder traversal sequence:   A, C, E, D, B, H, I, G, F
        Level-order traversal sequence: F, B, G, A, D, I, C, E, H
        """
        tree = trees.Tree(elem_id='tree1')
        tree.seed_node.elem_id = 'F'
        node_B = tree.seed_node.new_child('B')
        node_B.new_child('A')
        node_D = node_B.new_child('D')
        node_D.new_child('C')
        node_D.new_child('E')
        node_G = tree.seed_node.new_child('G')
        node_I = node_G.new_child('I')
        node_I.new_child('H')
        return tree

    def build_fixed_basic_tree2(self):
        """
        Returns a larger bifurcating basic tree:
        (   (   (   (   T1,
                        (   (   T2,
                                T3
                            )i5,
                            T4
                        )i4
                    )i3,
                    (   T5,
                        T6,
                        T13
                    )i6
                )i2,
                (   (   T7,
                        (   T8,
                            T9
                        )i9
                    )i8,
                    T10
                )i7
            )i1,
            T14,
            (   T11,
                T12
            )i10
        )i0;
        """
        tree = trees.Tree(elem_id='tree2')
        tree.seed_node.elem_id = 'i0'
        node_i1 = tree.seed_node.new_child('i1')
        node_i2 = node_i1.new_child('i2')
        node_i3 = node_i2.new_child('i3')
        node_i3.new_child('T1')
        node_i4 = node_i3.new_child('i4')
        node_i5 = node_i4.new_child('i5')
        node_i5.new_child('T2')
        node_i5.new_child('T3')
        node_i4.new_child('T4')
        node_i6 = node_i2.new_child('i6')
        node_i6.new_child('T5')
        node_i6.new_child('T6')
        node_i7 = node_i1.new_child('i7')
        node_i8 = node_i7.new_child('i8')
        node_i8.new_child('T7')
        node_i9 = node_i8.new_child('i9')
        node_i9.new_child('T8')
        node_i9.new_child('T9')
        node_i7.new_child('T10')
        tree.seed_node.new_child('T14')
        node_i10 = tree.seed_node.new_child('i10')
        node_i10.new_child('T11')
        node_i10.new_child('T12')
        node_i6.new_child('T13')
        return tree

    def setUp(self):        
        self.tree1 = self.build_fixed_basic_tree1()
        self.tree2 = self.build_fixed_basic_tree2()

    def show_edge(self, edge):
        show = []
        show.append('(')
        if edge.tail_node:
            show.append(edge.tail_node.elem_id)
        else:
            show.append('#')
        show.append(',')
        if edge.head_node:
            show.append(edge.head_node.elem_id)
        else:
            show.append('#')
        show.append(')')
        return ''.join(show)

    def test_leafiter(self):
        _LOG.debug('### LEAF ITERATION ###')
        result = [node.elem_id for node in self.tree1.leaf_iter()]
        expected = ["A", "C", "E", "H"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)
        result = [node.elem_id for node in self.tree2.leaf_iter()]
        expected = ["T1", "T2", "T3", "T4", "T5", "T6", "T13", "T7", "T8", "T9", "T10", "T14", "T11", "T12"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)

    def test_preorder(self):
        _LOG.debug('### PRE-ORDER ITERATION ###')    
        result = [node.elem_id for node in self.tree1.preorder_node_iter()]
        expected = ["F", "B", "A", "D", "C", "E", "G", "I", "H"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)
        result = [node.elem_id for node in self.tree2.preorder_node_iter()]
        expected = ["i0", "i1", "i2", "i3", "T1", "i4", "i5", "T2", "T3", "T4", "i6", "T5", "T6", "T13", "i7", "i8", "T7", "i9", "T8", "T9", "T10", "T14", "i10", "T11", "T12"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)

    def test_postorder(self):
        _LOG.debug('### POST-ORDER ITERATION ###')     
        result = [node.elem_id for node in self.tree1.postorder_node_iter()]
        expected = ["A", "C", "E", "D", "B", "H", "I", "G", "F"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)
        result = [node.elem_id for node in self.tree2.postorder_node_iter()]
        expected = ["T1", "T2", "T3", "i5", "T4", "i4", "i3", "T5", "T6", "T13", "i6", "i2", "T7", "T8", "T9", "i9", "i8", "T10", "i7", "i1", "T14", "T11", "T12", "i10", "i0"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)

    def test_levelorder(self):
        _LOG.debug('### LEVEL-ORDER ITERATION ###')     
        result = [node.elem_id for node in self.tree1.level_order_node_iter()]
        expected = ["F", "B", "G", "A", "D", "I", "C", "E", "H"]        
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)
        result = [node.elem_id for node in self.tree2.level_order_node_iter()]
        expected = ["i0", "i1", "T14", "i10", "i2", "i7", "T11", "T12", "i3", "i6", "i8", "T10", "T1", "i4", "T5", "T6", "T13", "T7", "i9", "i5", "T4", "T8", "T9", "T2", "T3"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        self.assertEqual(result, expected)
        
if __name__ == "__main__":
    unittest.main()

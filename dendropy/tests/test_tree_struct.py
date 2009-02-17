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
from dendropy.datasets import Dataset
import dendropy.tests
_LOG = get_logger("TreeStructureAndTraversal")

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
        tree = trees.Tree(oid='tree1')
        tree.seed_node.oid = 'F'
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
        tree = trees.Tree(oid='tree2')
        tree.seed_node.oid = 'i0'
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
            show.append(edge.tail_node.oid)
        else:
            show.append('#')
        show.append(',')
        if edge.head_node:
            show.append(edge.head_node.oid)
        else:
            show.append('#')
        show.append(')')
        return ''.join(show)

    def testLeafIter(self):
        _LOG.debug('### LEAF ITERATION ###')
        result = [node.oid for node in self.tree1.leaf_iter()]
        expected = ["A", "C", "E", "H"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected
        result = [node.oid for node in self.tree2.leaf_iter()]
        expected = ["T1", "T2", "T3", "T4", "T5", "T6", "T13", "T7", "T8", "T9", "T10", "T14", "T11", "T12"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected

    def testPreorder(self):
        _LOG.debug('### PRE-ORDER ITERATION ###')    
        result = [node.oid for node in self.tree1.preorder_node_iter()]
        expected = ["F", "B", "A", "D", "C", "E", "G", "I", "H"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected
        result = [node.oid for node in self.tree2.preorder_node_iter()]
        expected = ["i0", "i1", "i2", "i3", "T1", "i4", "i5", "T2", "T3", "T4", "i6", "T5", "T6", "T13", "i7", "i8", "T7", "i9", "T8", "T9", "T10", "T14", "i10", "T11", "T12"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected

    def testPostOrder(self):
        _LOG.debug('### POST-ORDER ITERATION ###')     
        result = [node.oid for node in self.tree1.postorder_node_iter()]
        expected = ["A", "C", "E", "D", "B", "H", "I", "G", "F"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected
        result = [node.oid for node in self.tree2.postorder_node_iter()]
        expected = ["T1", "T2", "T3", "i5", "T4", "i4", "i3", "T5", "T6", "T13", "i6", "i2", "T7", "T8", "T9", "i9", "i8", "T10", "i7", "i1", "T14", "T11", "T12", "i10", "i0"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected

    def testLevelOrder(self):
        _LOG.debug('### LEVEL-ORDER ITERATION ###')     
        result = [node.oid for node in self.tree1.level_order_node_iter()]
        expected = ["F", "B", "G", "A", "D", "I", "C", "E", "H"]        
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected
        result = [node.oid for node in self.tree2.level_order_node_iter()]
        expected = ["i0", "i1", "T14", "i10", "i2", "i7", "T11", "T12", "i3", "i6", "i8", "T10", "T1", "i4", "T5", "T6", "T13", "T7", "i9", "i5", "T4", "T8", "T9", "T2", "T3"]
        _LOG.debug("EXPECTED: %s" % (', '.join(result)))        
        _LOG.debug("  RESULT: %s" % (', '.join(expected)))        
        assert result == expected

    def testPostOrderAgain(self):
        n = '((((t49:0.0299,t41:0.0299):0.42017,(((t39:0.231767,(t7:0.213631,(t36:0.099739,((t12:0.015337,t26:0.015337):0.036523,t8:0.05186):0.047879):0.113892):0.018137):0.062275,((t43:0.240402,(((t28:0.157872,((t1:0.132296,t46:0.132296):0.016383,(t14:0.12343,t3:0.12343):0.025249):0.009193):0.01427,((t24:0.151816,(t25:0.106007,(t32:0.035143,(t47:0.025662,t6:0.025662):0.00948):0.070864):0.045809):0.013127,t40:0.164943):0.007199):0.010417,(t34:0.159553,(((t21:0.061001,t16:0.061001):0.001043,t44:0.062044):0.011229,(t18:0.001174,t5:0.001174):0.072098):0.086281):0.023006):0.057844):0.027253,(t27:0.066806,t50:0.066806):0.20085):0.026387):0.117517,(t23:0.259881,(t4:0.072245,((t45,t30):0.029839,t31:0.029839):0.042406):0.187635):0.151679):0.038511):0.069304,(t19:0.25004,t17:0.25004):0.269335):0.215312,((((t37:0.036909,t13:0.036909):0.244651,t2:0.28156):0.262824,((t15:0.11244,(t33:0.076665,t10:0.076665):0.035775):0.124232,t9:0.236671):0.307712):0.057112,(((t22:0.065248,t42:0.065248):0.037237,t38:0.102485):0.182409,(t48:0.187654,((t20:0.03615,(t29:0.01626,t11:0.01626):0.01989):0.039894,t35:0.076044):0.11161):0.09724):0.316601):0.133191);'
        d = Dataset()
        tree = d.trees_from_string(string=n, format="NEWICK")[0]
        tree.debug_check_tree()
        used = set()
        for node in tree.postorder_node_iter():
            ch = node.child_nodes()
            for c in ch:
                self.assertTrue(c in used)
            used.add(node)
        
if __name__ == "__main__":
    unittest.main()

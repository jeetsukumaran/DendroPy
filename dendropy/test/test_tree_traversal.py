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

class TestAgeOrderNodeIterator(unittest.TestCase):

    def setUp(self):
        self.tree = dendropy.Tree.get_from_string("""
               [&R] ((T1:7.83100438361,(((T2:0.765636721664,T3:0.765636721664)I05:1.5392268074,((T4:0.943318444552,T5:0.943318444552)I07:0.93553295287,T6:1.87885139742)I06:0.42601213164)I04:5.45683230098,((T7:2.49566105798,(T8:0.806690133601,(T9:0.52991878645,(T10:0.431923385995,(T11:0.125937061974,T12:0.125937061974)I13:0.305986324021)I12:0.097995400455)I11:0.27677134715)I10:1.68897092438)I09:2.18367702572,(T13:2.01769620826,(T14:1.60224165123,T15:1.60224165123)I15:0.415454557029)I14:2.66164187544)I08:3.08235774634)I03:0.0693085535678)I02:2.78993675296,((((T16:0.296004771084,T17:0.296004771084)I19:3.37487509214,T18:3.67087986323)I18:2.46214618289,((T19:4.53359864104,(((T20:1.53078041372,(T21:1.07630491497,T22:1.07630491497)I25:0.454475498748)I24:0.909377421616,(T23:2.11474343234,(T24:0.182766659126,T25:0.182766659126)I27:1.93197677321)I26:0.325414402998)I23:0.518984065864,(T26:1.16066496127,T27:1.16066496127)I28:1.79847693993)I22:1.57445673984)I21:1.36025142806,(T28:3.21016578878,T29:3.21016578878)I29:2.68368428032)I20:0.23917597702)I17:1.80259064407,T30:7.93561669019)I16:2.68532444638)I01:1.34961326574;""",
                "newick")
        for nd in self.tree:
            if nd.taxon is not None:
                nd.label = nd.taxon.label
            else:
                assert nd.label
        self.expected_leaves = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T11', 'T12',]
        self.expected_internal = ['I13', 'I27', 'I19', 'I12', 'I11', 'I05', 'I10', 'I07', 'I25', 'I28', 'I24', 'I15', 'I06', 'I14', 'I26', 'I04', 'I23', 'I09', 'I22', 'I29', 'I18', 'I21', 'I08', 'I20', 'I17', 'I03', 'I02', 'I16', 'I01']
        self.expected_all = self.expected_leaves + self.expected_internal

    def testAgeOrderIterAscendingIncludingLeaves(self):
        self.tree.calc_node_ages()
        results = [nd.label for nd in self.tree.age_order_node_iter(include_leaves=True)]
        self.assertEqual(set(self.expected_all), set(results))
        results = [i for i in results if i.startswith('I')]
        self.assertEqual(self.expected_internal, results)

    def testAgeOrderIterAscendingNotIncludingLeaves(self):
        self.tree.calc_node_ages()
        results = [nd.label for nd in self.tree.age_order_node_iter(include_leaves=False)]
        self.assertEqual(self.expected_internal, results)

    def testAgeOrderIterDescendingIncludingLeaves(self):
        self.tree.calc_node_ages()
        results = [nd.label for nd in self.tree.age_order_node_iter(include_leaves=True, descending=True)]
        self.assertEqual(set(self.expected_all), set(results))
        results = [i for i in results if i.startswith('I')]
        self.assertEqual(self.expected_internal[-1::-1], results)

    def testAgeOrderIterDescendingNotIncludingLeaves(self):
        self.tree.calc_node_ages()
        results = [nd.label for nd in self.tree.age_order_node_iter(include_leaves=False, descending=True)]
        self.assertEqual(self.expected_internal[-1::-1], results)

if __name__ == "__main__":
    unittest.main()

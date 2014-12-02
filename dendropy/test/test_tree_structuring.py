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
Tests native tree structuring routines.
"""

import unittest
from dendropy.test.support import pathmap
from dendropy.utility import messaging
from dendropy.test.support import extendedtest
from dendropy.test.support.datagen import RepeatedRandom
import dendropy
import re
import random

_LOG = messaging.get_logger(__name__)

class TestTreeLadderization(unittest.TestCase):

    def setUp(self):
        self.tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

    def clean_newick_str(self, s):
        """
        Strips out everything but the core tree statement characters from a
        NEWICK string.
        """
        return re.sub(r'[^()A-H,]', '', s)

    def testLadderizeLeft(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=True)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
                self.clean_newick_str("[&R] ((F,(G,H)),(A,(B,(C,(D,E)))));"))

    def testLadderizeRight(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=False)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
               self.clean_newick_str("[&R] (((((D,E),C),B),A),((G,H),F));"))

class TreeMidpointRootingTest(extendedtest.ExtendedTestCase):

    def testMidpointRooting(self):
        taxa = dendropy.TaxonSet()
        test_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.randomly-rooted.tre'),
                "nexus",
                taxon_set=taxa,
                as_rooted=True)
        expected_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.midpoint-rooted.tre'),
                "nexus",
                taxon_set=taxa,
                as_rooted=True)
        for idx, test_tree in enumerate(test_trees):
            expected_tree = expected_trees[idx]
            test_tree.reroot_at_midpoint(update_splits=True)
            self.assertEqual(test_tree.symmetric_difference(expected_tree), 0)
            for split in test_tree.split_edges:
                if test_tree.split_edges[split].head_node is test_tree.seed_node:
                    continue
                self.assertAlmostEqual(test_tree.split_edges[split].length, expected_tree.split_edges[split].length, 3)

class ResolvePolytomiesTestCase(unittest.TestCase):

    def verify_resolve_polytomies(self, tree_string, rng):
        if rng is None:
            return
        tree = dendropy.Tree.get_from_string(tree_string, "newick")
        for nd in tree:
            nd.edge.length = 100
        tree.resolve_polytomies(rng=rng)
        tree.encode_splits()
        tree.debug_check_tree()
        for nd in tree:
            if nd is tree.seed_node and not tree.is_rooted:
                self.assertEqual(len(nd._child_nodes), 3)
            elif len(nd._child_nodes) > 0:
                self.assertEqual(len(nd._child_nodes), 2)
        tree2 = dendropy.Tree.get_from_string(tree_string,
                "newick",
                taxon_set=tree.taxon_set)
        self.assertNotEqual(tree.symmetric_difference(tree2), 0)
        tree.collapse_unweighted_edges()
        tree.encode_splits()
        tree2.encode_splits()
        self.assertEqual(tree.symmetric_difference(tree2), 0)

    def test_resolve_polytomies_at_root(self):
        for tree_string in (
                "(a,b,c,d)e;",
                ):
            for rooting in ("[&R]", "[&U]"):
                tree_string2 = rooting + " " +  tree_string
                # cycle through rng period
                self.verify_resolve_polytomies(tree_string2, None)
                for x in range(1001):
                    rng = RepeatedRandom()
                    for i in range(x):
                        rng.uniform(0, 1)
                    self.verify_resolve_polytomies(tree_string2, rng)

    def xtest_resolve_polytomies(self):
        rng = random.Random()
        for tree_string in (
                "((((Homo:0.21,Bogus1:0.23,Pongo:0.21)N1:0.28,Bogus2:0.49,Macaca:0.49)N2:0.13,Bogus3:0.62,Ateles:0.62)N3:0.38,Galago:1.00)N4:0.0;",
                "(((t52,t62,(t2,t58,(t32,(t55,t28,t39,t17,t4,t44,t25)internal6,t26,t9,t48,(t41,t45)internal7)internal5,t56)internal4,t54,((t18,t14,t34)internal9,(t49,t22,t50,t27,t16,t40,t6,t19)internal10,t13,(t51,t35,t61,t53,t43)internal11)internal8,((t42,t5,t7,t33,t30,t21,t47,t38)internal13,t23,t1,t11,t46)internal12,(t63,t3,t37,t59)internal14)internal3,t57,t64,t31,(t12,(t60,t24,t10,t15,t20)internal16,t36)internal15,t8)internal2,t29)internal1;",
                "((t13,t37,t19,((t21,t44)internal4,(t61,t46,t4,t8,t63)internal5,t23,t28)internal3,t52,t64,(t39,t40,t24)internal6,(t54,t62,t15,t55,t51)internal7)internal2,(t3,(t53,t33,(t47,t9,t25)internal10,t45,t18,t27,t17)internal9,(t10,t22)internal11,(t59,t20,t12,t57,t56,t38,t7,t11)internal12,((t1,t31,t43,t36,t34,t14,t49)internal14,t2,t41,t50,((t30,t32,t58,t6)internal16,t60,t16)internal15,t26,t48,t35)internal13,t42,t29,t5)internal8)internal1;",
                "((t7,t13,(t42,t51,t20,t26,(t21,t18,t16)internal4)internal3,(t6,t48,(t23,t33,t34,t15,t2,t25)internal6,(t64,t45,t49,t3,t55,t31,(t19,t47,t38,t35)internal8,t14)internal7,t9,(t5,(t62,t50)internal10,t54,t32,(t40,t8,t58,t60,t10,t30)internal11)internal9)internal5)internal2,t28,(t4,t57,t52,t43,t46,t39)internal12,t63,t11,(t44,(t12,t22,t36,t29,t24,t1,t17,t56)internal14,(t41,t59,t53,t61)internal15,t37)internal13,t27)internal1;",
                "(t10,t42,(t54,t12,t40,(t30,(t55,t3,t22,t56)internal4,t39)internal3,t51,((t41,(t28,t52,t24,t14,t49,t38,t36,(t35,t34,t13,t9,t59,t58)internal8)internal7,t5,t45,t17,(t23,(t11,t53,t57,t19,t26)internal10)internal9,t63,(t62,t29,t18,t20,t27,t43)internal11)internal6,t4,(t25,t7,(((t47,t61,t21,t64)internal15,t44,t46,t15,t6,t37,t48)internal14,t50,t60)internal13)internal12)internal5,t31)internal2,t1,t32,t33,(t2,t8)internal16,t16)internal1;",
                "(t60,(((t29,t35,t5,t59,t4,t9,(t6,t25,t37,t44,t54)internal5)internal4,(t53,t23,t28,t48,((t16,t46,t26,t10,t52)internal8,t30,t17,t51,t40)internal7,t31,t3)internal6,t14,t7,t2,t22,(t45,t56,t20,t36,t43,t47)internal9,t8)internal3,t18,(t39,(t61,t27,t21,t58,t24)internal11,t50,(t64,t38,(t32,t11,t49,t63)internal13)internal12,t13,((t57,t34,t15,t55,t19)internal15,t1,t42)internal14,t33)internal10,(t62,t12,t41)internal16)internal2)internal1;",
                "((t12,t16,(t51,t24,t27,(t38,t44,t52,t6,t9,t53,t20,t18)internal4,(t31,t47,t56,t60)internal5)internal3,t34,t42)internal2,t5,t62,((t61,t64,t59,t15,t48,t2,t35)internal7,t57,(t14,t28,t40,t22,t58,t7,(t25,t4)internal9,((t8,t29,t21)internal11,t17,t1,t46,(t33,t50,t11)internal12,t13,t45,t37)internal10)internal8,(t26,t10,(t3,t30,t32,(t19,t36,t39,t43)internal15,t23,t54,t55)internal14)internal13,t41,t49,t63)internal6)internal1;",
                "(((t31,t19,t12,(t55,t25)internal4,t11,t60,(t46,t8,t56)internal5)internal3,(t52,t18)internal6,(((t24,t16,t7)internal9,t57,t40,(t2,t9,t50,t37,t43,t20,t15,t22)internal10,t42,(t47,t28,t58,t10)internal11)internal8,t21,(t14,(t51,t53,t26,t35,t49)internal13,t33,t62,t34,t54,t17)internal12,t5,t61,t4)internal7,t3,t36,t45,t63)internal2,(t27,t32,t30,t48)internal14,(t13,t64,(t23,t59,t41,t38,t44,t1)internal16,t29,t39,t6)internal15)internal1;",
                "(((t8,t40,t25,t36,t37,t11)internal3,t29,t31,t51,t3,((t61,t43,t63,t50,t23,t52,t24,t30)internal5,t39,t44,t58,((t15,t64,t9,t28,(t5,t34,t38,t22,(t33,t12,t35,t42,t10,t2,t27,t45)internal9)internal8,t7,t21,t60)internal7,t41,(t48,t17,t14)internal10)internal6,t13,t18)internal4)internal2,((t47,t16,t49,(t53,t26)internal13,(t19,t20,(t55,t54)internal15,t32,t4)internal14,t59,t57,t6)internal12,t46,t56,t1)internal11,t62)internal1;",
                "(t33,t39,t18,t35,((t19,t62,t55,(t41,t14,(t1,t36,t16,t38,t25,t59,t34)internal5)internal4,t61,t50)internal3,(((t2,t48,t22)internal8,t28,t37,t47,(t60,t30)internal9,t27)internal7,t12,(t31,t21,(t3,t5,t45,t53)internal11,(t23,t54,t20,t4,t64,t56,(t58,t13,t26,t11,t57,t44,t42,t46)internal13)internal12,t40)internal10)internal6,t7,(t52,(t43,t24)internal15,(t49,t29,t63,t32)internal16,t9,t6,t8,t15)internal14,t51)internal2,t10,t17)internal1;",
                "(t35,((t60,t47,t58,t26,t9)internal3,t3,((t13,t25)internal5,t1,((t6,t32,(t53,t7,t64,t34,t18,t23,t30,t33)internal8,t55,(t48,t20,t12,t4,t38,t28,t36)internal9)internal7,t8,((t57,t40,t52,t31,t43,t39)internal11,t59,(t37,t16,t27,(t44,t41,t10)internal13,t50)internal12,t24,t63)internal10,(t5,t56,t61,t29)internal14,t21,t49)internal6)internal4,((t17,t42,t62,t15)internal16,(t19,t2,t51,(t22,t14,(t45,t54,t46)internal19)internal18)internal17,t11)internal15)internal2)internal1;",
                ):
            for rooting in ("[&R]", "[&U]"):
                tree_string2 = rooting + " " +  tree_string
                self.verify_resolve_polytomies(tree_string2, rng)

if __name__ == "__main__":
    unittest.main()


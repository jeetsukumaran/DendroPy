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
Splits encoding tests.
"""

import sys
import os
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.utility import messaging
import dendropy
from dendropy import treesplit
from dendropy import treecalc
from dendropy import treemanip

_LOG = messaging.get_logger(__name__)

class CladeMaskTest(unittest.TestCase):

    def runTest(self):
        tree_list = dendropy.TreeList(
            stream=StringIO("""((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);"""),
            schema="newick")
        for i in tree_list:
            _LOG.debug(i.get_indented_form())
            treesplit.encode_splits(i)
            _LOG.debug(i.get_indented_form(splits=True))
            i.debug_check_tree(splits=True, logger_obj=_LOG)
        root1 = tree_list[0].seed_node
        root1e = root1.edge
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask), range(6))
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, one_based=True), range(1,7))
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, mask=21, one_based=True), [1, 3, 5])
        self.assertEqual(treesplit.split_to_list(root1e.split_bitmask, mask=21), [0, 2, 4])
        self.assertEqual(treesplit.count_bits(root1e.split_bitmask), 6)

        fc1 = root1.child_nodes()[0]
        fc1e = fc1.edge
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask), [0, 1])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, one_based=True), [1, 2])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, mask=0x15, one_based=True), [1])
        self.assertEqual(treesplit.split_to_list(fc1e.split_bitmask, mask=0x15), [0])
        self.assertEqual(treesplit.count_bits(fc1e.split_bitmask), 2)

class CountBitsTest(unittest.TestCase):

    def runTest(self):
        self.assertEqual(treesplit.count_bits(21), 3)

class LowestBitTest(unittest.TestCase):

    def runTest(self):
        for n, expected in enumerate([0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16]):
            self.assertEqual(treesplit.lowest_bit_only(n), expected)

class IsTrivialTest(unittest.TestCase):

    def runTest(self):
        y = True
        n = False
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0xF))
        for i, r in enumerate([y, y, y, n, y, n, n, n, y, n, n, n, n, n, n, y, y, n, n, n, n, n, n, y, n, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0x1F))
                              #0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, treesplit.is_trivial_split(i, 0x17))

if __name__ == "__main__":
    unittest.main()

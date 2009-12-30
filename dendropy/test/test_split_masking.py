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

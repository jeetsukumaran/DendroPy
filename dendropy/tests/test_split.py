#! /usr/bin/env python

############################################################################
##  test_tree_gens.py
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
Tests tree distances.
"""

import random
import unittest
import math
from dendropy import get_logger
import dendropy.tests
_LOG = get_logger("TreeGenerationAndSimulation")

from dendropy import dataio
from dendropy.tests.util_for_testing import assert_approx_equal, assert_vec_approx_equal, assert_mat_approx_equal

### MODULE THAT WE ARE TESTING ###
from dendropy.splits import *
### MODULE THAT WE ARE TESTING ###

class SplitTest(unittest.TestCase):

    def testCladeMasks(self):
        dataset = dataio.trees_from_newick([
                                       "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       ])
        tree_list = [i[0] for i in dataset.trees_blocks]
        for i in tree_list:
            _LOG.debug(i.get_indented_form())
            encode_splits(i)
            _LOG.debug(i.get_indented_form(splits=True))
            i.debug_check_tree(splits=True, logger_obj=_LOG)
        root1 = tree_list[0].seed_node
        root1e = root1.edge
        self.assertEqual(split_to_list(root1e.clade_mask), range(6))
        self.assertEqual(split_to_list(root1e.clade_mask, one_based=True), range(1,7))
        self.assertEqual(split_to_list(root1e.clade_mask, mask=21, one_based=True), [1, 3, 5])
        self.assertEqual(split_to_list(root1e.clade_mask, mask=21), [0, 2, 4])
        self.assertEqual(count_bits(root1e.clade_mask), 6)

        fc1 = root1.child_nodes()[0]
        fc1e = fc1.edge
        self.assertEqual(split_to_list(fc1e.clade_mask), [0, 1])
        self.assertEqual(split_to_list(fc1e.clade_mask, one_based=True), [1, 2])
        self.assertEqual(split_to_list(fc1e.clade_mask, mask=0x15, one_based=True), [1])
        self.assertEqual(split_to_list(fc1e.clade_mask, mask=0x15), [0])
        self.assertEqual(count_bits(fc1e.clade_mask), 2)


    def testCountBits(self):
        self.assertEqual(count_bits(21), 3)
    
    def testLowestBitOnly(self):
        for n, expected in enumerate([0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16]):
            self.assertEqual(lowest_bit_only(n), expected)

    def testIsTrivial(self):
        y = True
        n = False
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, is_trivial_split(i, 0xF))
        for i, r in enumerate([y, y, y, n, y, n, n, n, y, n, n, n, n, n, n, y, y, n, n, n, n, n, n, y, n, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, is_trivial_split(i, 0x1F))
                              #0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1 
        for i, r in enumerate([y, y, y, n, y, n, n, y, y, y, y, n, y, n, n, y, y, n, n, y, n, y, y, y, y, n, n, y, n, y, y, y, ]):
            self.assertEqual(r, is_trivial_split(i, 0x17))
if __name__ == "__main__":
    unittest.main()

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

    def testEuclideanDist(self):
        dataset = dataio.trees_from_newick([
                                       "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:0.161175,t6:0.161175):0.392293,((t2:0.075411,(t4:0.104381,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);",
                                       ])
        tree_list = [i[0] for i in dataset.trees_blocks]
        #print "\n".join([str(i) for i in tree_list])
        for i in tree_list:
            encode_splits(i, taxa_block= dataset.taxa_blocks[0])
        root1 = tree_list[0].seed_node
        root1e = root1.edge
        fc1 = root1.child_nodes()[0]
        fc1e = fc1.edge
        self.assertEqual(split_to_list(root1e.split_mask), range(6))
        self.assertEqual(split_to_list(root1e.split_mask, one_based=True), range(1,7))
        self.assertEqual(split_to_list(root1e.split_mask, mask=21, one_based=True), [1, 3, 5])
        self.assertEqual(split_to_list(root1e.split_mask, mask=21), [0, 2, 4])
        self.assertEqual(count_bits(root1e.split_mask), 6)

        self.assertEqual(split_to_list(fc1e.split_mask), range(2))
        self.assertEqual(split_to_list(fc1e.split_mask, one_based=True), range(1,3))
        self.assertEqual(split_to_list(fc1e.split_mask, mask=21, one_based=True), [1])
        self.assertEqual(split_to_list(fc1e.split_mask, mask=21), [0])
        self.assertEqual(count_bits(fc1e.split_mask), 2)


    def testCountBits(self):
        self.assertEqual(count_bits(21), 3)

if __name__ == "__main__":
    unittest.main()

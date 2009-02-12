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
_LOG = get_logger("tree diste   ")

from dendropy import dataio
from dendropy.splits import encode_splits
from dendropy.tests.util_for_testing import assert_approx_equal, assert_vec_approx_equal, assert_mat_approx_equal

### MODULE THAT WE ARE TESTING ###
from dendropy import treedists
### MODULE THAT WE ARE TESTING ###

class TreeDistTest(unittest.TestCase):

    def testEuclideanDist(self):
        d = dataio.trees_from_newick([
                                       "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:0.161175,t6:0.161175):0.392293,((t2:0.075411,(t4:0.104381,t1:0.075411):1):0.065840,t3:0.170221):0.383247);",
                                       "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);",
                                       ])
        tree_list = [i[0] for i in d.trees_blocks]
        #print "\n".join([str(i) for i in tree_list])
        for i in tree_list:
            encode_splits(i)
        assert_approx_equal(treedists.euclidean_distance(tree_list[0], tree_list[1]), 2.0)
        assert_approx_equal(treedists.euclidean_distance(tree_list[0], tree_list[2]), math.sqrt(2.0))
        assert_approx_equal(treedists.euclidean_distance(tree_list[0], tree_list[3]), 0.97103099999999998)

        assert_approx_equal(treedists.euclidean_distance(tree_list[1], tree_list[2]), math.sqrt(6.0))
        assert_approx_equal(treedists.euclidean_distance(tree_list[1], tree_list[3]), 2.2232636377544162)

        assert_approx_equal(treedists.euclidean_distance(tree_list[2], tree_list[3]), 1.000419513484718)
        
    def testSymmDiff(self):
        newick = "((t5,t6),((t4,(t2,t1)),t3));"
        d = dataio.trees_from_newick([newick])
        ref = d.trees_blocks[0][0]
        taxa_block = d.taxa_blocks[0]


        encode_splits(ref)

        o_newick = "((t1,t2),((t4,(t5,t6)),t3));"
        o_tree = dataio.trees_from_newick([o_newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(o_tree)

        self.assertEqual(treedists.symmetric_difference(o_tree, ref), 2)


if __name__ == "__main__":
    unittest.main()

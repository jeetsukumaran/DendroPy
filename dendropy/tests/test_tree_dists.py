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
from dendropy import get_logger
import dendropy.tests
_LOG = get_logger("TreeGenerationAndSimulation")

from dendropy import dataio
from dendropy.splits import encode_splits
from dendropy.tests.util_for_testing import assert_approx_equal, assert_vec_approx_equal, assert_mat_approx_equal

### MODULE THAT WE ARE TESTING ###
from dendropy import treedists
### MODULE THAT WE ARE TESTING ###

class TreeDistTest(unittest.TestCase):

    def testEuclideanDist(self):
        return
        newick1 = "(t1:1,t2:2,(t3:3,t4:4):5)"  #"((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        newick3 = "(t1:1,t3:3,(t2:2,t4:4):5)"  # "((t5:0.161175,t6:0.161175):0.392293,((t2:0.075411,(t4:0.104381,t1:0.075411):1):0.065840,t3:0.170221):0.383247);"
        tree1 = dataio.trees_from_string(string=newick1, format="NEWICK")[0]
        encode_splits(tree1)
        tree3 = dataio.trees_from_string(string=newick3, format="NEWICK", taxa_block=tree1)[0]
        encode_splits(tree3)
        d13 = treedists.euclidean_distance(tree1, tree3)
        assert_approx_equal(d13, 7.0710678118654755)
        pass


        newick2 = "((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        newick4 = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        tree2 = dataio.trees_from_string(string=newick2, format="NEWICK")[0]
        encode_splits(tree2)
        tree4 = dataio.trees_from_string(string=newick4, format="NEWICK")[0]
        encode_splits(tree4)
        return
        d12 = treedists.euclidean_distance(tree1, tree2)
        assert_approx_equal(d12, 2.0)
        d14 = treedists.euclidean_distance(tree1, tree4)
        assert_approx_equal(d14, 2.0)
        d23 = treedists.euclidean_distance(tree2, tree3)
        assert_approx_equal(d23, 2.0)
        d24 = treedists.euclidean_distance(tree2, tree4)
        assert_approx_equal(d24, 2.0)
        d34 = treedists.euclidean_distance(tree3, tree4)
        assert_approx_equal(d34, 2.0)
        


if __name__ == "__main__":
    unittest.main()

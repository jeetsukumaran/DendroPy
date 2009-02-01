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
from dendropy.treemanip import *
### MODULE THAT WE ARE TESTING ###

class TreeManipTest(unittest.TestCase):

    def testCollapseEdge(self):
        tree = dataio.trees_from_newick(["((t5,t6),((t4,(t2,t1)),t3));"]).trees_blocks[0][0]
        root = tree.seed_node
        fc = root.child_nodes()[0]
        collapse_edge(fc.edge)
        self.assertEqual(str(tree), "(t5,t6,((t4,(t2,t1)),t3))")

    def testCollapseClade(self):
        tree = dataio.trees_from_newick(["((t5,t6),((t4,(t2,t1)),t3));"]).trees_blocks[0][0]
        root = tree.seed_node
        root_children = root.child_nodes()
        fc = root_children[0]
        collapse_clade(fc)
        self.assertEqual(str(tree), "((t5,t6),((t4,(t2,t1)),t3))")
        fc1 = root_children[1]
        fc1children = fc1.child_nodes()
        fc1child = fc1children[0]
        collapse_clade(fc1child)
        self.assertEqual(str(tree), "((t5,t6),((t4,t2,t1),t3))")
        collapse_clade(fc1)
        self.assertEqual(str(tree), "((t5,t6),(t4,t2,t1,t3))")
        collapse_clade(root)
        self.assertEqual(str(tree), "(t5,t6,t4,t2,t1,t3)")
        tree = dataio.trees_from_newick(["((t5,t6),((t4,(t2,t1)),t3));"]).trees_blocks[0][0]
        root = tree.seed_node
        collapse_clade(root)
        self.assertEqual(str(tree), "(t5,t6,t4,t2,t1,t3)")
        

if __name__ == "__main__":
    unittest.main()

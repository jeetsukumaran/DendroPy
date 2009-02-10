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
import copy
from dendropy import get_logger
from dendropy.treedists import symmetric_difference
from dendropy.splits import encode_splits
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

    def testReroot(self):
        newick = "((t5,t6),((t4,(t2,t1)),t3));"
        d = dataio.trees_from_newick([newick])
        tree = d.trees_blocks[0][0]
        taxa_block = d.taxa_blocks[0]
        ref = dataio.trees_from_newick([newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(ref, taxa_block=taxa_block)

        o_newick = "((t2, t1),((t4,(t5,t6)),t3));"
        o_tree = dataio.trees_from_newick([o_newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(o_tree, taxa_block=taxa_block)
        self.assertEqual(symmetric_difference(o_tree, ref), 2)

        taxa_labels = ["t%d" % i for i in xrange(1,7)]
        for leaf_name in taxa_labels:
            f = lambda x : x.label == leaf_name
            nd = tree.find_taxon_node(f)
            tree.to_outgroup_position(nd)
            r_newick =  str(tree)
            r_tree = dataio.trees_from_newick([r_newick], taxa_block=taxa_block).trees_blocks[0][0]
            encode_splits(r_tree, taxa_block=taxa_block)
            self.assertEqual(symmetric_difference(r_tree, ref), 0)

if __name__ == "__main__":
    unittest.main()

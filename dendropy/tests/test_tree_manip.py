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
from dendropy.splits import encode_splits, split_to_list, count_bits, lowest_bit_only
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
        self.assertEqual(str(tree), "((t5,t6),((t4,(t2,t1)),t3))")
        fc = root.child_nodes()[0]
        collapse_edge(fc.edge)
        tree.debug_check_tree(logger_obj=_LOG)
        self.assertEqual(str(tree), "(t5,t6,((t4,(t2,t1)),t3))")

    def testCollapseClade(self):
        tree = dataio.trees_from_newick(["(t5,t6,((t4,(t2,t1)),t3));"]).trees_blocks[0][0]
        encode_splits(tree)
        root = tree.seed_node
        root_children = root.child_nodes()
        fc = root_children[0]
        collapse_clade(fc)
        tree.debug_check_tree(splits=True)
        self.assertEqual(str(tree), "(t5,t6,((t4,(t2,t1)),t3))")
        fc2 = root_children[2]
        fc2children = fc2.child_nodes()
        t124child = fc2children[0]
        collapse_clade(t124child)
        tree.debug_check_tree(logger_obj=_LOG)
        self.assertEqual(str(tree), "(t5,t6,((t4,t2,t1),t3))")
        collapse_clade(fc2)
        tree.debug_check_tree(logger_obj=_LOG)
        self.assertEqual(str(tree), "(t5,t6,(t4,t2,t1,t3))")
        collapse_clade(root)
        tree.debug_check_tree(logger_obj=_LOG)
        tree.debug_check_tree(logger_obj=_LOG)
        self.assertEqual(str(tree), "(t5,t6,t4,t2,t1,t3)")
        tree = dataio.trees_from_newick(["((t5,t6),((t4,(t2,t1)),t3));"]).trees_blocks[0][0]
        root = tree.seed_node
        collapse_clade(root)
        tree.debug_check_tree(logger_obj=_LOG)
        self.assertEqual(str(tree), "(t5,t6,t4,t2,t1,t3)")

    def testReroot(self):
        newick = "((t5,t6),((t4,(t2,t1)),t3));"
        d = dataio.trees_from_newick([newick])
        tree = d.trees_blocks[0][0]
        taxa_block = d.taxa_blocks[0]
        ref = dataio.trees_from_newick([newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(ref)

        o_newick = "((t2, t1),((t4,(t5,t6)),t3));"
        o_tree = dataio.trees_from_newick([o_newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(o_tree)
        self.assertEqual(symmetric_difference(o_tree, ref), 2)

        taxa_labels = ["t%d" % i for i in xrange(1,7)]
        for leaf_name in taxa_labels:
            f = lambda x : x.label == leaf_name
            nd = tree.find_taxon_node(f)
            tree.to_outgroup_position(nd)
            r_newick =  str(tree)
            r_tree = dataio.trees_from_newick([r_newick], taxa_block=taxa_block).trees_blocks[0][0]
            encode_splits(r_tree)
            self.assertEqual(symmetric_difference(r_tree, ref), 0)

    def testRerootSplits(self):
        newick = "((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));"
        d = dataio.trees_from_newick([newick])
        tree = d.trees_blocks[0][0]
        taxa_block = d.taxa_blocks[0]
        ref = dataio.trees_from_newick([newick], taxa_block=taxa_block).trees_blocks[0][0]
        encode_splits(tree)
        encode_splits(ref)
        r = tree.seed_node
        curr_n = r.child_nodes()[1]
        former_mask = curr_n.edge.clade_mask
        tm = r.edge.clade_mask
        nbits = count_bits(tm)
        from dendropy.splits import split_as_string
        
        tree.reroot_at(curr_n, splits=True, delete_deg_two=False)
        
        new_root = tree.seed_node
        self.assertEqual(tm, new_root.edge.clade_mask)
        self.assertEqual(True, new_root is curr_n)
        self.assertEqual(True, r.parent_node is curr_n)
        flipped = (~(r.edge.clade_mask)) & tm
        self.assertEqual(True, (former_mask == r.edge.clade_mask) or (flipped == former_mask))

if __name__ == "__main__":
    unittest.main()

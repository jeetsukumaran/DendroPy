
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
Tests of tree metrics.
"""

import random
import unittest
import math
from cStringIO import StringIO

import dendropy
from dendropy.treecalc import fitch_down_pass
from dendropy.test.support import pathmap

from dendropy import treecalc

class FitchTest(unittest.TestCase):

    def testPScore(self):
        expected_scores = [370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 370, 671, 670, 678, 687, 633, 675, 689, 668, 652, 644]
        dataset = dendropy.DataSet(stream=open(pathmap.char_source_path("apternodus.chars.nexus"), "rU"),
                                   format='NEXUS')
        dataset.read(stream=open(pathmap.tree_source_path("apternodus.tre"), "rU"), 
                     format='NEXUS',
                     taxon_set=dataset.taxon_sets[0])
        char_mat = dataset.char_matrices[0]
        taxa_to_state_set_map = char_mat.create_taxon_to_state_set_map()
        tree_list = dataset.tree_lists[0]
        self.assertEqual(len(expected_scores), len(tree_list))
        for n, tree in enumerate(tree_list):
            node_list = tree.postorder_node_iter()
            pscore = fitch_down_pass(node_list, taxa_to_state_set_map=taxa_to_state_set_map)
            self.assertEqual(expected_scores[n], pscore)

if __name__ == "__main__":
    unittest.main()


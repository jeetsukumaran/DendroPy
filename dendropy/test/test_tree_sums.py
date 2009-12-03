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
Tests of summarization.
"""

import unittest
import dendropy
from dendropy import treecalc
from dendropy import treesum
from dendropy import treesplit
from dendropy.test.support import pathmap

class TestConsensusTree(unittest.TestCase):

    def setUp(self):
#        x = dendropy.Tree.get_from_path(
#                pathmap.tree_source_path("pythonidae_cytb.mb.con"),
#                format="nexus",
#                index=0)
#        self.tree_list = dendropy.TreeList(taxon_set=x.taxon_set)
        self.tree_list = dendropy.TreeList()
        for t in xrange(1, 5):
            tf = pathmap.tree_source_path('pythonidae_cytb.mb.run%d.t' % t)
            self.tree_list.read_from_path(tf, 'nexus', tree_offset=25)
        self.mb_con_tree = dendropy.Tree.get_from_path(
                pathmap.tree_source_path("pythonidae_cytb.mb.con"),
                format="nexus",
                index=0,
                taxon_set=self.tree_list.taxon_set)
        self.mb_con_tree.encode_splits()

    def testConsensus(self):
        con_tree = self.tree_list.consensus(min_freq=0.50, trees_splits_encoded=False, support_label_decimals=2)
        con_tree.encode_splits()

#        print self.tree_list.taxon_set.description(2)
#        print
#        print self.mb_con_tree.taxon_set.description(2)
#        print
#        print con_tree.taxon_set.description(2)

#        print "#NEXUS"
#        print "Begin trees;"
#        print "tree 1 = ", self.mb_con_tree.as_string("newick"), ";"
#        print "tree 2 = ", con_tree.as_string("newick"),";"
#        print "end;"
        self.assertEqual(treecalc.symmetric_difference(self.mb_con_tree, con_tree), 0)
        self.assertEqual(len(con_tree.split_edges), len(self.mb_con_tree.split_edges))
        sd = self.tree_list.split_distribution
        for split in self.mb_con_tree.split_edges:
            edge1 = self.mb_con_tree.split_edges[split]
            edge2 = con_tree.split_edges[split]
#            print con_tree.taxon_set.split_bitmask_string(split), edge1.head_node.label, edge2.head_node.label, sd.split_frequencies[split]
            if edge1.head_node.label and edge2.head_node.label:
                s1 = float(edge1.head_node.label)
                s2 = round(float(edge2.head_node.label), 2)
                self.assertAlmostEqual(s1, s2, 2)

if __name__ == "__main__":
    unittest.main()

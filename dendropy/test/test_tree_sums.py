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
from dendropy import treesum
from dendropy import treesplit
from dendropy.test.support import pathmap

class TestConsensusTree(unittest.TestCase):

    def setUp(self):
        self.tree_list = dendropy.TreeList()
        for t in xrange(1, 4):
            tf = pathmap.tree_source_path('pythonidae_cytb.mb.run%d.t' % t)
            self.tree_list.read_from_path(tf, 'nexus')
        self.mb_con_tree = dendropy.Tree.get_from_path(
                pathmap.tree_source_path("pythonidae_cytb.mb.con"),
                format="nexus",
                from_index=0,
                taxon_set=self.tree_list.taxon_set)

    def testConsensus(self):
        con_tree = self.tree_list.consensus(min_freq=0.5, trees_splits_encoded=False)
        con_tree.encode_splits()

if __name__ == "__main__":
    unittest.main()

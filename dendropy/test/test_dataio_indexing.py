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
Tests of data indexing.
"""

import unittest
import dendropy

class TestTreeIndexing(unittest.TestCase):

    def setUp(self):
        self.trees1 = """
#NEXUS

BEGIN TREES
    TITLE = TREES1;
    TREE 0_0 = ((A,B),(C,D));
    TREE 0_1 = ((A,B),(C,D));
    TREE 0_2 = ((A,B),(C,D));
END;

BEGIN TREES;
    TITLE = TREES2;
    TREE 1_0 = ((A,B),(C,D));
    TREE 1_1 = ((A,B),(C,D));
    TREE 1_2 = ((A,B),(C,D));
END;

BEGIN TREES;
    TITLE = TREES3;
    TREE 2_0 = ((A,B),(C,D));
    TREE 2_1 = ((A,B),(C,D));
    TREE 2_2 = ((A,B),(C,D));
END;

BEGIN TREES
    TITLE = TREES4;
    TREE 3_0 = ((A,B),(C,D));
    TREE 3_1 = ((A,B),(C,D));
    TREE 3_2 = ((A,B),(C,D));
END;

"""

    def testDefaultTreeIndexing(self):
        t = dendropy.Tree.get_from_string(self.trees1, "nexus")
        self.assertEqual(t.label, '0 0')

    def testUnifiedTreeIndexing(self):
        t = dendropy.Tree.get_from_string(self.trees1, "nexus", tree_offset=7)
        self.assertEqual(t.label, '2 1')

    def testCollectionTreeIndexing(self):
        t = dendropy.Tree.get_from_string(self.trees1, "nexus", collection_offset=3, tree_offset=2)
        self.assertEqual(t.label, '3 2')

    def testDefaultTreeOutOfRange(self):
        self.assertRaises(IndexError, dendropy.Tree.get_from_string, self.trees1, "nexus", tree_offset=99)

    def testCollectionTreeOutOfRange(self):
        self.assertRaises(IndexError, dendropy.Tree.get_from_string, self.trees1, "nexus", collection_offset=7, tree_offset=0)



if __name__ == "__main__":
    unittest.main()

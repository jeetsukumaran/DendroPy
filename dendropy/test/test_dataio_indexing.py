#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests of data indexing.
"""

import unittest
import dendropy

class IndexingTestCase(unittest.TestCase):

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

class TestTreeIndexing(IndexingTestCase):

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

class TestTreeListIndexing(IndexingTestCase):

    def testDefaultCollection(self):
        t = dendropy.TreeList.get_from_string(self.trees1, "nexus")
        self.assertEqual(len(t), 12)

    def testDefaultCollectionTreeOffset(self):
        t = dendropy.TreeList.get_from_string(self.trees1, "nexus", tree_offset=7)
        self.assertEqual(len(t), 5)
        self.assertEqual(t[0].label, '2 1')

    def testTreeCollectionIndexing(self):
        t = dendropy.TreeList.get_from_string(self.trees1, "nexus", collection_offset=2)
        self.assertEqual(len(t), 3)
        self.assertEqual(t[0].label, '2 0')

    def testTreeCollectionIndexingTreeIndexing(self):
        t = dendropy.TreeList.get_from_string(self.trees1, "nexus", collection_offset=2, tree_offset=1)
        self.assertEqual(len(t), 2)
        self.assertEqual(t[0].label, '2 1')

    def testDefaultCollectionTreeOutOfRange(self):
        self.assertRaises(IndexError, dendropy.TreeList.get_from_string, self.trees1, "nexus", tree_offset=99)

    def testCollectionOutOfRange(self):
        self.assertRaises(IndexError, dendropy.TreeList.get_from_string, self.trees1, "nexus", collection_offset=7, tree_offset=0)

    def testCollectionTreeOutOfRange(self):
        self.assertRaises(IndexError, dendropy.TreeList.get_from_string, self.trees1, "nexus", collection_offset=2, tree_offset=4)

if __name__ == "__main__":
    unittest.main()

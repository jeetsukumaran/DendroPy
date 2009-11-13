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
Tests creation, reading, update, deletion of TreeList objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import datatest
from dendropy.test.support import datagen
from dendropy.test.support import pathmap
import dendropy

class TreeListCreateTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.tree_list1 = datagen.reference_tree_list()
        self.tree_list1_stream = StringIO(self.tree_list1.as_string("newick"))

    def testDeepCopyTreeListFromTreeListSameTaxa(self):
        tree_list2 = dendropy.TreeList(self.tree_list1, taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True)

    def testDeepCopyTreeListFromTreeListSameTaxa(self):
        tree_list2 = dendropy.TreeList([dendropy.Tree(t) for t in self.tree_list1], taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True)

    def testShallowCopyTreeListFromListSameTaxa(self):
        src_trees = [t for t in self.tree_list1]
        tree_list2 = dendropy.TreeList(src_trees, taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=False)

    def testDeepCopyTreeListFromTreeListDifferentTaxa(self):
        tree_list2 = dendropy.TreeList(self.tree_list1)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True)

    def testDeepCopyTreeListFromTreeListDifferentTaxa(self):
        tree_list2 = dendropy.TreeList([dendropy.Tree(t) for t in self.tree_list1])
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True)

    def testTooManyPosArgs(self):
        self.assertRaises(error.TooManyArgumentsError, dendropy.TreeList, self.tree_list1, dendropy.TreeList())

    def testMultipleSources(self):
        self.assertRaises(error.MultipleInitializationSourceError, dendropy.TreeList, self.tree_list1, stream=self.tree_list1_stream, format="newick")

    def testTreeListFromFileSameTaxa(self):
        tree_list2 = dendropy.TreeList(stream=self.tree_list1_stream, format="newick", taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True)

    def testTreeListFromFileDifferentTaxa(self):
        tree_list2 = dendropy.TreeList([dendropy.Tree(t) for t in self.tree_list1])
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True)

    def testTreeListFromFileNoFormatSpecification(self):
        self.assertRaises(error.UnspecifiedFormatError, dendropy.TreeList, stream=self.tree_list1_stream)

    def testTreeListFromFileNoKeywords(self):
        self.assertRaises(ValueError, dendropy.TreeList, self.tree_list1_stream)

    def testFromFileFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        s = pathmap.tree_source_path('reference.trees.nexus')
        tree_list2 = dendropy.TreeList.from_file(open(s, "rU"), "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

    def testFromPathFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        s = pathmap.tree_source_path('reference.trees.nexus')
        tree_list2 = dendropy.TreeList.from_path(s, "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

    def testFromStringFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        tree_list2 = dendropy.TreeList.from_string(tree_list1.as_string('nexus'), "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

if __name__ == "__main__":
    unittest.main()

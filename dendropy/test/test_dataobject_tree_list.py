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
        self.tree_list1_stream = StringIO(self.tree_list1.as_string("nexus"))

    def testDeepCopyTreeListFromTreeListSameTaxa(self):
        tree_list2 = dendropy.TreeList(self.tree_list1, taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True, ignore_label=True)

    def testDeepCopyTreeListFromTreeListSameTaxa(self):
        tree_list2 = dendropy.TreeList([dendropy.Tree(t) for t in self.tree_list1], taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True, ignore_label=True)

    def testShallowCopyTreeListFromListSameTaxa(self):
        src_trees = [t for t in self.tree_list1]
        tree_list2 = dendropy.TreeList(src_trees, taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=False, ignore_label=True)

    def testDeepCopyTreeListFromTreeListDifferentTaxa(self):
        tree_list2 = dendropy.TreeList(self.tree_list1)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True, ignore_label=True)

    def testDeepCopyTreeListFromTreeListDifferentTaxa(self):
        tree_list2 = dendropy.TreeList([dendropy.Tree(t) for t in self.tree_list1])
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True, ignore_taxon_order=True)

    def testTooManyPosArgs(self):
        self.assertRaises(error.TooManyArgumentsError, dendropy.TreeList, self.tree_list1, dendropy.TreeList())

    def testMultipleSources(self):
        self.assertRaises(error.MultipleInitializationSourceError, dendropy.TreeList, self.tree_list1, stream=self.tree_list1_stream, schema="newick")

    def testTreeListFromFileSameTaxa(self):
        tree_list2 = dendropy.TreeList(stream=self.tree_list1_stream, schema="nexus", taxon_set=self.tree_list1.taxon_set)
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=False, equal_oids=False, distinct_trees=True)

    def testTreeListFromFileDifferentTaxa(self):
        tree_list2 = dendropy.TreeList(stream=self.tree_list1_stream, schema="nexus")
        self.assertDistinctButEqual(self.tree_list1, tree_list2, distinct_taxa=True, equal_oids=False, distinct_trees=True)

    def testTreeListFromFileNoFormatSpecification(self):
        self.assertRaises(error.UnspecifiedSchemaError, dendropy.TreeList, stream=self.tree_list1_stream)

    def testTreeListFromFileNoKeywords(self):
        self.assertRaises(ValueError, dendropy.TreeList, self.tree_list1_stream)

    def testFromFileFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree_list2 = dendropy.TreeList.get_from_stream(open(s, "rU"), "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

    def testFromPathFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree_list2 = dendropy.TreeList.get_from_path(s, "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

    def testFromStringFactoryDistinctTaxa(self):
        tree_list1 = datagen.reference_tree_list()
        tree_list2 = dendropy.TreeList.get_from_string(tree_list1.as_string('nexus'), "nexus")
        self.assertDistinctButEqual(tree_list1, tree_list2, distinct_taxa=True)

if __name__ == "__main__":
    unittest.main()

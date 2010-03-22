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
Tests creation, reading, update, deletion of Tree objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import datatest
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
import dendropy

class TreeCreateTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.tree1 = datagen.four_taxon_tree1()
        self.tree1_newick_str = self.tree1.as_newick_string(include_internal_labels=True)

    def testTreeFromTreeSameTaxa(self):
        tree2 = dendropy.Tree(self.tree1)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromTreeSetOidAndLabelSameTaxa(self):
        tree2 = dendropy.Tree(self.tree1, oid="TREE2", label="TREE2")
        self.assertEqual(tree2.oid, "TREE2")
        self.assertEqual(tree2.label, "TREE2")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromTreeDistinctTaxa(self):
        taxa = dendropy.TaxonSet()
        tree2 = dendropy.Tree(self.tree1, taxon_set=taxa)
        self.assertIs(tree2.taxon_set, taxa)
        self.assertIsNot(tree2.taxon_set, self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromTreeWithExtraPosArgs(self):
        self.assertRaises(TypeError, dendropy.Tree, self.tree1, "dummy")

    def testTreeFromInvalidObjectPosArgs(self):
        self.assertRaises(error.InvalidArgumentValueError, dendropy.Tree, object())

    def testTreeFromInvalidIterablePosArgs(self):
        self.assertRaises(error.InvalidArgumentValueError, dendropy.Tree, "abcde")

    def testTreeFromFileTooManyPosArgs(self):
        self.assertRaises(error.TooManyArgumentsError, dendropy.Tree, StringIO(self.tree1_newick_str), "newick")

    def testTreeFromFileKeywordArgsDistinctTaxa(self):
        tree2 = dendropy.Tree(stream=StringIO(self.tree1_newick_str), schema="newick")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromFileKeywordArgsSameTaxa(self):
        tree2 = dendropy.Tree(stream=StringIO(self.tree1_newick_str), schema="newick", taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromFileMixedPosAndKeywordArgs(self):
        self.assertRaises(error.MultipleInitializationSourceError, dendropy.Tree, self.tree1, tream=StringIO(self.tree1_newick_str), schema="newick")

    def testTreeFromTreeWithExtraKeywordArgsOK(self):
        tree2 = dendropy.Tree(self.tree1, stream=None, schema=None)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromFilePosArgsWithNoFormat(self):
        self.assertRaises(error.UnspecifiedSchemaError, dendropy.Tree, stream=StringIO(self.tree1_newick_str), taxon_set=self.tree1.taxon_set)

    def testTreeFromMultipleSource(self):
        self.assertRaises(error.MultipleInitializationSourceError,
                dendropy.Tree, \
                StringIO(self.tree1_newick_str),
                stream=StringIO(self.tree1_newick_str),
                schema="newick",
                taxon_set=self.tree1.taxon_set)

    def testTreeFromReadDistinctTaxa(self):
        tree2 = dendropy.Tree()
        tree2.read_from_string(self.tree1_newick_str, "newick")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromReadSameTaxa(self):
        tree2 = dendropy.Tree()
        tree2.read_from_string(self.tree1_newick_str, "newick", taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromNewickReadIndexed(self):
        nstr = "(A,(B,(C,D))); ((A,C),(B,D)); %s; (A,(C,(B,D))); ((A,D),(B,C));" % self.tree1_newick_str
        tree2 = dendropy.Tree()
        tree2.read_from_string(nstr, "newick", tree_offset=2, taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromNewickFileIndexed(self):
        nstr = "(A,(B,(C,D))); ((A,C),(B,D)); %s; (A,(C,(B,D))); ((A,D),(B,C));" % self.tree1_newick_str
        tree2 = dendropy.Tree(stream=StringIO(nstr), schema="newick", tree_offset=2, taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testFromFileFactoryDistinctTaxa(self):
        tree_list = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree = dendropy.Tree.get_from_stream(open(s, "rU"), "nexus", tree_offset=2)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=True)

    def testFromPathFactoryDistinctTaxa(self):
        tree_list = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree = dendropy.Tree.get_from_path(s, "nexus", tree_offset=2)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=True)

    def testFromStringFactoryDistinctTaxa(self):
        tree_list = datagen.reference_tree_list()
        tree = dendropy.Tree.get_from_string(tree_list.as_string('nexus'), "nexus", tree_offset=2)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=True)

    def testFromFileFactorySameTaxa(self):
        tree_list = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree = dendropy.Tree.get_from_stream(open(s, "rU"), "nexus", tree_offset=2, taxon_set=tree_list.taxon_set)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=False)

    def testFromPathFactorySameTaxa(self):
        tree_list = datagen.reference_tree_list()
        s = pathmap.tree_source_path(datagen.reference_trees_filename(schema="nexus"))
        tree = dendropy.Tree.get_from_path(s, "nexus", tree_offset=2, taxon_set=tree_list.taxon_set)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=False)

    def testFromStringFactorySameTaxa(self):
        tree_list = datagen.reference_tree_list()
        tree = dendropy.Tree.get_from_string(tree_list.as_string('nexus'), "nexus", tree_offset=2, taxon_set=tree_list.taxon_set)
        self.assertDistinctButEqual(tree_list[2], tree, distinct_taxa=False)


if __name__ == "__main__":
    unittest.main()

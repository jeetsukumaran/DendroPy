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
Tests creation, reading, update, deletion of Tree and TreeList objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import errors
from dendropy.test.support import framework
from dendropy.test.support import datagen
import dendropy

class TreeInstantiationTest(framework.DataObjectVerificationTestCase):

    def setUp(self):
        self.tree1 = datagen.get_standard_four_taxon_tree()
        self.tree1_newick_str = self.tree1.to_newick_str(include_internal_labels=True)

    def testTreeFromStandard(self):
        node_oids = [nd.oid for nd in self.tree1.postorder_node_iter()]
        self.assertEqual(node_oids, ['a', 'b', 'i1', 'c', 'd', 'i2', 'root'])
        tax_labels = [nd.taxon.label for nd in self.tree1.postorder_node_iter() if nd.taxon is not None]
        self.assertEqual(tax_labels, ['A', 'B', 'C', 'D'])

    def testTreeFromTreeSameTaxa(self):
        tree2 = dendropy.Tree(self.tree1)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromTreeSetOidAndLabelSameTaxa(self):
        tree2 = dendropy.Tree(self.tree1, oid="TREE2", label="TREE2")
        self.assertEqual(tree2.oid, "TREE2")
        self.assertEqual(tree2.label, "TREE2")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromTreeWithExtraKeywordArgs(self):
        self.assertRaises(TypeError, dendropy.Tree, self.tree1, format="newick")

    def testTreeFromTreeWithExtraPosArgs(self):
        self.assertRaises(TypeError, dendropy.Tree, self.tree1, "dummy")

    def testTreeFromInvalidObjectPosArgs(self):
        self.assertRaises(TypeError, dendropy.Tree, object())

    def testTreeFromInvalidObjectPosArgsWithKeywords(self):
        self.assertRaises(TypeError, dendropy.Tree, "xx", stream=StringIO(self.tree1_newick_str), format="newick")

    def testTreeFromFileKeywordArgsDistinctTaxa(self):
        tree2 = dendropy.Tree(istream=StringIO(self.tree1_newick_str), format="newick")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromFileKeywordArgsSameTaxa(self):
        tree2 = dendropy.Tree(istream=StringIO(self.tree1_newick_str), format="newick", taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromFilePosArgsDistinctTaxa(self):
        tree2 = dendropy.Tree(StringIO(self.tree1_newick_str), "newick")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromFilePosArgsSameTaxa(self):
        tree2 = dendropy.Tree(StringIO(self.tree1_newick_str), "newick", taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromFileMixedPosAndKeywordArgsDistinctTaxa(self):
        tree2 = dendropy.Tree(StringIO(self.tree1_newick_str), format="newick")
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=True, equal_oids=False)

    def testTreeFromFileMixedPosAndKeywordArgsSameTaxa(self):
        tree2 = dendropy.Tree(StringIO(self.tree1_newick_str), format="newick", taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromFilePosArgsWithNoFormat(self):
        self.assertRaises(errors.UnspecifiedFormatError, dendropy.Tree, istream=StringIO(self.tree1_newick_str), taxon_set=self.tree1.taxon_set)

    def testTreeFromMultipleSource(self):
        self.assertRaises(TypeError,
                dendropy.Tree, \
                StringIO(self.tree1_newick_str),
                istream=StringIO(self.tree1_newick_str),
                format="newick",
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
        tree2.read_from_string(nstr, "newick", index=2, taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

    def testTreeFromNewickFileIndexed(self):
        nstr = "(A,(B,(C,D))); ((A,C),(B,D)); %s; (A,(C,(B,D))); ((A,D),(B,C));" % self.tree1_newick_str
        tree2 = dendropy.Tree(StringIO(nstr), "newick", index=2, taxon_set=self.tree1.taxon_set)
        self.assertDistinctButEqual(self.tree1, tree2, distinct_taxa=False, equal_oids=False)

#    def test_treelist_init_from_newick(self):
#
#        newick_str = "((A,B),(C,D)); ((A,C),(B,D)); (A,(B,(C,D))); (A,(C,(B,D)));"
#
#        # from file, using keywords
#        tl1 = dendropy.TreeList(istream=StringIO(newick_str), format="newick", oid="t1")
#        self.assertTrue(tl1.oid == "t1", "'%s'" % tl1.oid)
#        self.assertEqual(len(tl1), 4)
#        self.assertEqual(len(tl1.taxon_set), 4)
#        for t in tl1:
#            t.debug_check_tree(_LOG)
#            self.assertTrue(tl1.taxon_set is t.taxon_set)
#
#        # test copying
#        tl2 = dendropy.TreeList(tl1)
#        self.assertTrue(tl2.taxon_set is tl1.taxon_set)
#        self.assertTrue(tl2.oid != tl1.oid)
#        self.assertTrue(tl2.label == tl1.label)
#        self.assertEqual(len(tl2.taxon_set), 4)
#        self.assertTrue(tl2.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
#        self.assertEqual(len(tl1), len(tl2))
#        for ti, t1 in enumerate(tl1):
#            t2 = tl2[ti]
#            self.assertTrue(t2 is t1)
#
#        # from file, args
#        tl3 = dendropy.TreeList(StringIO(newick_str), "newick", taxon_set=tl1.taxon_set)
#        self.assertTrue(tl3.taxon_set is tl1.taxon_set)
#
#        # from file, mixed
#        tl4 = dendropy.TreeList(StringIO(newick_str), format="newick", taxon_set=tl1.taxon_set)
#        self.assertTrue(tl4.taxon_set is tl1.taxon_set)
#
#        # read from string
#        tl5 = dendropy.TreeList()
#        tl5.read_from_string(newick_str, format="newick")

if __name__ == "__main__":
    unittest.main()

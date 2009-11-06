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
from dendropy.utility import messaging
from dendropy.test.support import framework
from dendropy.test.support import datagen
_LOG = messaging.get_logger(__name__)

import dendropy

class TreeInstantiationTest(framework.DataObjectVerificationTestCase):

    def testTreeFromStandard(self):
        tree = datagen.get_standard_four_taxon_tree()
        node_oids = [nd.oid for nd in tree.postorder_node_iter()]
        self.assertEqual(node_oids, ['a', 'b', 'i1', 'c', 'd', 'i2', 'root'])
        tax_labels = [nd.taxon.label for nd in tree.postorder_node_iter() if nd.taxon is not None]
        self.assertEqual(tax_labels, ['A', 'B', 'C', 'D'])

    def testTreeFromTree(self):
        tree1 = datagen.get_standard_four_taxon_tree()
        tree2 = dendropy.Tree(tree1)
        self.assertDistinctButEqualTrees(tree1, tree2, distinct_taxa=False, equal_oids=False)

#        datacheck.compare_individual_trees(tree1, tree2, tester=self, distinct_taxa=True, distinct_oids=False)

#    def test_tree_init_from_newick(self):
#        newick_str = "((A,B),(C,D));"
#
#        # from file, using keywords
#        t1 = dendropy.Tree(istream=StringIO(newick_str), format="newick", oid="t1")
#        self.assertTrue(t1.oid == "t1", "'%s'" % t1.oid)
#        t1.debug_check_tree(_LOG)
#
#        # test copying
#        t2 = dendropy.Tree(t1)
#        t2.debug_check_tree(_LOG)
#        self.compare_tree_copies(t1, t2)
#
#        # from file, args
#        t3 = dendropy.Tree(StringIO(newick_str), "newick", taxon_set=t1.taxon_set)
#        t3.debug_check_tree(_LOG)
#        self.assertTrue(t3.taxon_set is t1.taxon_set)
#
#        # from file, mixed
#        t4 = dendropy.Tree(StringIO(newick_str), format="newick", taxon_set=t1.taxon_set)
#        t4.debug_check_tree(_LOG)
#        self.assertTrue(t4.taxon_set is t1.taxon_set)
#
#        # read from string
#        t5 = dendropy.Tree()
#        t5.read_from_string(newick_str, format="newick")
#
#    def test_tree_init_from_nexus(self):
#
#        nexus_str = """\
#NEXUS
#begin taxa;
#    dimensions ntax=4;
#    taxlabels
#        A
#        B
#        C
#        D
#    ;
#end;
#begin trees;
#    translate
#        1 A,
#        2 B,
#        3 C,
#        4 D;
#    tree 1 = ((A,B)i1, (C,D)i2)root;
#end;
#"""
#        # NEXUS
#        t1 = dendropy.Tree(istream=StringIO(nexus_str), format="nexus")
#        t2 = dendropy.Tree(StringIO(nexus_str), format="nexus")
#        t3 = dendropy.Tree(StringIO(nexus_str), "nexus")
#        t4 = dendropy.Tree()
#        t4.read_from_string(nexus_str, "nexus")
#        t5 = dendropy.Tree()
#        t5.read_from_file(StringIO(nexus_str), "nexus")
#        t6 = dendropy.Tree(StringIO(nexus_str), "nexus", taxon_set=t1.taxon_set)
#        self.assertTrue(t6.taxon_set is t1.taxon_set)
#        t7 = dendropy.Tree()
#        t7.read_from_string(nexus_str, "nexus", taxon_set=t1.taxon_set)
#        self.assertTrue(t7.taxon_set is t1.taxon_set)
#        for tx in (t1, t2, t3, t4, t5):
#            tx.debug_check_tree(_LOG)
#            tx.taxon_set.is_mutable = False
#            self.assertEqual(len(tx.taxon_set), 4, str([t.label for t in tx.taxon_set]))
#            self.assertTrue(tx.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
#        t8 = dendropy.Tree(t1)
#        self.compare_tree_copies(t1, t8)
#
#    def compare_tree_copies(self, t1, t2):
#        self.assertTrue(t2 is not t1)
#        self.assertTrue(t2.taxon_set is t1.taxon_set)
#        self.assertTrue(t2.seed_node is not t1.seed_node)
#        self.assertNotEqual(t1.oid, t2.oid)
#        self.assertEqual(t1.label, t2.label)
#        self.assertNotEqual(t1.oid, t2.oid)
#        t1_nodes = [nd for nd in t1.postorder_node_iter()]
#        t2_nodes = [nd for nd in t2.postorder_node_iter()]
#        for ndi, nd1 in enumerate(t1_nodes):
#            nd2 = t2_nodes[ndi]
#            self.assertTrue(nd1 is not nd2)
#            self.assertNotEqual(nd1.oid, nd2.oid)
#            self.assertEqual(nd1.label, nd2.label)
#            self.assertTrue(nd1.taxon is nd2.taxon, "%s vs. %s" % (repr(nd1.taxon), repr(nd2.taxon)))
#        t1_edges = [e for e in t1.postorder_edge_iter()]
#        t2_edges = [e for e in t2.postorder_edge_iter()]
#        for ei, e1 in enumerate(t1_edges):
#            e2 = t2_edges[ei]
#            self.assertTrue(e1 is not e2)
#            self.assertNotEqual(e1.oid, e2.oid)
#            self.assertEqual(e1.label, e2.label)
#
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

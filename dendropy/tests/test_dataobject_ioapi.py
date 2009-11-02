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
Tests composition and traversal of trees.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import messaging
import dendropy.tests
from dendropy import Tree, TreeList, TaxonSet

_LOG = messaging.get_logger(__name__)

class TreeInstantiationTest(unittest.TestCase):

#     def test_tree_init_from_newick(self):
#
#         newick_str = "((A,B),(C,D));"
#
#         # from file, using keywords
#         t1 = Tree(istream=StringIO(newick_str), format="newick", oid="t1")
#         self.assertTrue(t1.oid == "t1", "'%s'" % t1.oid)
#         t1.debug_check_tree(_LOG)
#
#         # test copying
#         t2 = Tree(t1)
#         t2.debug_check_tree(_LOG)
#         self.compare_trees(t1, t2)
#
#         # from file, args
#         t3 = Tree(StringIO(newick_str), "newick", taxon_set=t1.taxon_set)
#         t3.debug_check_tree(_LOG)
#         self.assertTrue(t3.taxon_set is t1.taxon_set)
#
#         # from file, mixed
#         t4 = Tree(StringIO(newick_str), format="newick", taxon_set=t1.taxon_set)
#         t4.debug_check_tree(_LOG)
#         self.assertTrue(t4.taxon_set is t1.taxon_set)
#
#         # read from string
#         t5 = Tree()
#         t5.read_from_string(newick_str, format="newick")
#
#     def test_tree_init_from_nexus(self):
#
#         nexus_str = """\
# #NEXUS
# begin taxa;
#     dimensions ntax=4;
#     taxlabels
#         A
#         B
#         C
#         D
#     ;
# end;
# begin trees;
#     translate
#         1 A,
#         2 B,
#         3 C,
#         4 D;
#     tree 1 = ((A,B)i1, (C,D)i2)root;
# end;
# """
#         # NEXUS
#         t1 = Tree(istream=StringIO(nexus_str), format="nexus")
#         t2 = Tree(StringIO(nexus_str), format="nexus")
#         t3 = Tree(StringIO(nexus_str), "nexus")
#         t4 = Tree()
#         t4.read_from_string(nexus_str, "nexus")
#         t5 = Tree()
#         t5.read_from_file(StringIO(nexus_str), "nexus")
#         t6 = Tree(StringIO(nexus_str), "nexus", taxon_set=t1.taxon_set)
#         self.assertTrue(t6.taxon_set is t1.taxon_set)
#         t7 = Tree()
#         t7.read_from_string(nexus_str, "nexus", taxon_set=t1.taxon_set)
#         self.assertTrue(t7.taxon_set is t1.taxon_set)
#         for tx in (t1, t2, t3, t4, t5):
#             tx.debug_check_tree(_LOG)
#             tx.taxon_set.is_mutable = False
#             self.assertEqual(len(tx.taxon_set), 4, str([t.label for t in tx.taxon_set]))
#             self.assertTrue(tx.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
#         t8 = Tree(t1)
#         self.compare_trees(t1, t8)
#
#     def compare_trees(self, t1, t2):
#         self.assertTrue(t2 is not t1)
#         self.assertTrue(t2.taxon_set is t1.taxon_set)
#         self.assertTrue(t2.seed_node is not t1.seed_node)
#         self.assertNotEqual(t1.oid, t2.oid)
#         self.assertEqual(t1.label, t2.label)
#         self.assertNotEqual(t1.oid, t2.oid)
#         t1_nodes = [nd for nd in t1.postorder_node_iter()]
#         t2_nodes = [nd for nd in t2.postorder_node_iter()]
#         for ndi, nd1 in enumerate(t1_nodes):
#             nd2 = t2_nodes[ndi]
#             self.assertTrue(nd1 is not nd2)
#             self.assertNotEqual(nd1.oid, nd2.oid)
#             self.assertEqual(nd1.label, nd2.label)
#             self.assertTrue(nd1.taxon is nd2.taxon, "%s vs. %s" % (repr(nd1.taxon), repr(nd2.taxon)))
#         t1_edges = [e for e in t1.postorder_edge_iter()]
#         t2_edges = [e for e in t2.postorder_edge_iter()]
#         for ei, e1 in enumerate(t1_edges):
#             e2 = t2_edges[ei]
#             self.assertTrue(e1 is not e2)
#             self.assertNotEqual(e1.oid, e2.oid)
#             self.assertEqual(e1.label, e2.label)

    def test_treelist_init_from_newick(self):

        newick_str = "((A,B),(C,D)); ((A,C),(B,D)); (A,(B,(C,D))); (A,(C,(B,D)));"

        # from file, using keywords
        tl1 = TreeList(istream=StringIO(newick_str), format="newick", oid="t1")
        self.assertTrue(tl1.oid == "t1", "'%s'" % tl1.oid)
        self.assertEqual(len(tl1), 4)
        self.assertEqual(len(tl1.taxon_set), 4)
        for t in tl1:
            t.debug_check_tree(_LOG)
            self.assertTrue(tl1.taxon_set is t.taxon_set)

        # test copying
        tl2 = TreeList(tl1)
        self.assertTrue(tl2.taxon_set is tl1.taxon_set)
        self.assertTrue(tl2.oid != tl1.oid)
        self.assertTrue(tl2.label == tl1.label)
        self.assertEqual(len(tl2.taxon_set), 4)
        self.assertTrue(tl2.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
        for ti, t1 in enumerate(tl1):
            t2 = tl2[ti]
            self.assertTrue(t2 is t1)

        # from file, args
        tl3 = TreeList(StringIO(newick_str), "newick", taxon_set=tl1.taxon_set)
        self.assertTrue(tl3.taxon_set is tl1.taxon_set)

        # from file, mixed
        tl4 = TreeList(StringIO(newick_str), format="newick", taxon_set=tl1.taxon_set)
        self.assertTrue(tl4.taxon_set is tl1.taxon_set)

        # read from string
        tl5 = TreeList()
        tl5.read_from_string(newick_str, format="newick")

if __name__ == "__main__":
    unittest.main()

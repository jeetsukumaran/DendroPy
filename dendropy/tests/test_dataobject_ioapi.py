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
from dendropy import Tree, TaxonSet

_LOG = messaging.get_logger(__name__)

class TreeInstantiationTest(unittest.TestCase):

    def test_tree_init(self):

        newick_str = "((A,B),(C,D));"
        nexus_str = """\
#NEXUS
begin taxa;
    dimensions ntax=4;
    taxlabels =
        A
        B
        C
        D
    ;
end;
begin trees;
    translate
        1 A,
        2 B,
        3 C,
        4 D;
    tree 1 = ((A,B)i1, (C,D)i2)root;
end;
"""

        # from file, using keywords
        t1 = Tree(istream=StringIO(newick_str), format="newick", oid="t1")
        self.assertTrue(t1.oid == "t1", "'%s'" % t1.oid)
        t1.debug_check_tree(_LOG)

        # test copying
        t2 = Tree(t1)
        t2.debug_check_tree(_LOG)
        self.compare_trees(t1, t2)

        # from file, args
        t3 = Tree(StringIO(newick_str), "newick")
        t3.debug_check_tree(_LOG)

        # from file, mixed
        t4 = Tree(StringIO(newick_str), format="newick")
        t4.debug_check_tree(_LOG)

        # NEXUS
        t5 = Tree(istream=StringIO(nexus_str), format="nexus")
        t6 = Tree(StringIO(nexus_str), format="nexus")
        t7 = Tree(StringIO(nexus_str), "nexus")
        for tx in (t5, t6, t7):
            tx.is_mutable = False
            self.assertEqual(len(tx.taxon_set), str([t.label for t in tx.taxon_set]))
            self.assertTrue(tx.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
            self.assertTrue(tx.taxon_set.has_taxon(label="A"))
            self.assertTrue(tx.taxon_set.has_taxa(taxa=self.taxon_set))
        t8 = Tree(t5)

    def compare_trees(self, t1, t2):
        self.assertTrue(t2 is not t1)
        self.assertTrue(t2.taxon_set is t1.taxon_set)
        self.assertTrue(t2.seed_node is not t1.seed_node)
        self.assertTrue(t1.oid != t2.oid)
        self.assertTrue(t1.label == t2.label)
        t1.oid = "t1"
        t2.oid = "t2"
        self.assertTrue(t1.oid != t2.oid)
        self.assertTrue(t1.oid == "t1", "'%s'" % t1.oid)
        self.assertTrue(t2.oid == "t2")
        if t1.label is None:
            self.assertTrue(t1.label is None)
        else:
            self.assertTrue(t1.label is not t2.label)
            self.assertTrue(t1.label == t2.label)
        t1.label = "t1"
        t2.label = "t2"
        self.assertTrue(t1.label != t2.label)
        self.assertTrue(t1.label == "t1", "'%s'" % t1.label)
        self.assertTrue(t2.label == "t2")
        t1_nodes = [nd for nd in t1.postorder_node_iter()]
        t2_nodes = [nd for nd in t2.postorder_node_iter()]
        for ndi, nd1 in enumerate(t1_nodes):
            nd2 = t2_nodes[ndi]
            self.assertTrue(nd1 is not nd2)
            self.assertTrue(nd1.oid != nd2.oid)
            self.assertTrue(nd1.label == nd2.label)
            if nd1.label is None:
                self.assertTrue(nd1.label is None)
            else:
                self.assertTrue(nd1.label is not nd2.label)
            self.assertTrue(nd1.taxon is nd2.taxon, "%s vs. %s" % (repr(nd1.taxon), repr(nd2.taxon)))
        t1_edges = [e for e in t1.postorder_edge_iter()]
        t2_edges = [e for e in t2.postorder_edge_iter()]
        for ei, e1 in enumerate(t1_edges):
            e2 = t2_nodes[ei]
            self.assertTrue(e1 is not e2)
            self.assertTrue(e1.oid != e2.oid)
            self.assertTrue(e1.label == e2.label)
            if e1.label is None:
                self.assertTrue(e1.label is None)
            else:
                self.assertTrue(e1.label is not e2.label)


if __name__ == "__main__":
    unittest.main()

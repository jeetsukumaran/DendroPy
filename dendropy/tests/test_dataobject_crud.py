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
Tests create, read, update, and delete methods of phylogenetic data objects.
"""

import random
import unittest
from cStringIO import StringIO
from dendropy.utility import messaging
from dendropy.tests import datacompare
import dendropy.tests
import dendropy

_LOG = messaging.get_logger(__name__)

class DatasetInstantiationTest(unittest.TestCase):

    def test_dataset_nexus_dna(self):
        nexus_str = """\
#NEXUS

Begin Taxa;
  dimensions ntax=4;
  taxlabels A B C D;
End;

Begin Characters;
   dimensions nchar=8;
   format datatype=dna missing=? gap=-;
   matrix
   A    AAAAAAAA
   B    ACTGACTG
   C    CCCCCCCC
   D    GGGGGGGG
   ;
End;

Begin Trees;
    tree 1 = [&R] ((A,B),(C,D));
    tree 2 = [&R] ((A,C),(B,D));
End;
"""
        d1 = dendropy.Dataset(StringIO(nexus_str), "nexus")
        self.assertEqual(len(d1.taxon_sets), 1)
        self.assertEqual(len(d1.taxon_sets[0]), 4)
        self.assertEqual(['A','B','C','D'], [t.label for t in d1.taxon_sets[0]])
        self.assertEqual(len(d1.char_arrays), 1)
        self.assertEqual(len(d1.tree_lists), 1)

        d2 = dendropy.Dataset(d1)

        self.assertTrue(d2.taxon_sets is not d1.taxon_sets)
        self.assertEqual(len(d2.taxon_sets), len(d1.taxon_sets))
        for tsi, ts1 in enumerate(d1.taxon_sets):
            ts2 = d2.taxon_sets[tsi]
            self.assertTrue(ts2 is not ts1)
            self.assertEqual(len(ts1), len(ts2))
            for i, t1 in enumerate(ts1):
                t2 = ts2[i]
                self.assertTrue(t2 is not t1)
                self.assertEqual(t2.label, t1.label)
                self.assertEqual(t2.oid, t1.oid)

        self.assertTrue(d2.tree_lists is not d1.tree_lists)
        self.assertEqual(len(d2.tree_lists), len(d1.tree_lists))
        for tli, tl1 in enumerate(d1.tree_lists):
            tl2 = d2.tree_lists[tli]
            self.assertTrue(tl2 is not tl1)
            self.assertEqual(len(tl2), len(tl1))
            self.assertTrue(tl2.taxon_set is not tl1.taxon_set)
            self.assertTrue(tl2.taxon_set in d2.taxon_sets)
            self.assertTrue
            for ti, t1 in enumerate(tl1):
                t2 = tl2[ti]
                self.assertTrue(t2 is not t1)
                self.assertEqual(t2.oid, t1.oid)
                self.assertEqual(t2.label, t1.label)
                self.assertTrue(t2.taxon_set is not t1.taxon_set)
                self.assertTrue(t2.taxon_set is tl2.taxon_set)


# class CharArrayInstantiationTest(unittest.TestCase):
#
#     def test_dna(self):
#
#         taxon_set = dendropy.TaxonSet(['A','B','C','D'])
#         ca1 = dendropy.DnaCharacterArray(taxon_set=taxon_set)
#         self.assertTrue(ca1.taxon_set is taxon_set)
#
#         def get_dna_cells():
#             cells = [dendropy.CharacterDataCell(value=s) for s in dendropy.DNA_STATE_ALPHABET.get_states(symbols="ACGTACGT")]
#             return cells
#
#         dna_col = dendropy.ColumnType(state_alphabet=dendropy.DNA_STATE_ALPHABET)
#         ca1.column_types.append(dna_col)
#         ca1[taxon_set[0]] = dendropy.CharacterDataVector(self.get_dna_cells())
#         ca1[taxon_set[1]] = dendropy.CharacterDataVector(self.get_dna_cells())
#         ca1[taxon_set[2]] = dendropy.CharacterDataVector(self.get_dna_cells())
#         ca1[taxon_set[3]] = dendropy.CharacterDataVector(self.get_dna_cells())
#
#         for t, v in ca1.items():
#             for c in v:
#                 c.column_type = dna_col
#                 self.assertTrue(c.value in dendropy.DNA_STATE_ALPHABET, c.value)
#
#         ca2 = dendropy.DnaCharacterArray(ca1)
#         self.assertTrue(ca1.taxon_set is ca2.taxon_set)
#         self.assertTrue(ca1.state_alphabets is ca2.state_alphabets)
#         self.assertEqual(ca1.state_alphabets, ca2.state_alphabets)
#         self.assertTrue(ca1.default_state_alphabet is ca2.default_state_alphabet)
#         self.assertEqual(len(ca2.column_types), 1)
#         self.assertTrue(ca2.column_types[0] is not ca1.column_types[0])
#         self.assertTrue(ca2.column_types[0].state_alphabet is ca1.column_types[0].state_alphabet)
#         for t, v1 in ca1.items():
#             v2 = ca2[t]
#             self.assertTrue(v1 is not v2)
#             for i, c1 in enumerate(v1):
#                 c2 = v2[i]
#                 self.assertTrue(c1 is not c2)
#                 self.assertTrue(c1.column_type is not c2.column_type)
#                 self.assertTrue(c1.column_type.state_alphabet is c2.column_type.state_alphabet)
#                 self.assertTrue(c1.value is c2.value, [id(c1.value), id(c2.value)])
#                 self.assertTrue(c1.value in dendropy.DNA_STATE_ALPHABET)
#                 self.assertTrue(c2.value in dendropy.DNA_STATE_ALPHABET)
#
#     def get_dna_cells(self):
#         cells = [dendropy.CharacterDataCell(value=s) for s in dendropy.DNA_STATE_ALPHABET.get_states(symbols="ACGTACGT")]
#         return cells
#
#     def test_standard(self):
#
#         taxon_set = dendropy.TaxonSet(['A','B','C','D'])
#         ca1 = dendropy.StandardCharacterArray(taxon_set=taxon_set)
#         self.assertTrue(ca1.taxon_set is taxon_set)
#
#         sa1 = self.get_standard_state_alphabet("012")
#         sa2 = self.get_standard_state_alphabet("XYZ")
#         sa3 = self.get_standard_state_alphabet("JKL")
#         ca1.state_alphabets = [sa1, sa2, sa3]
#         col_012 = dendropy.ColumnType(state_alphabet=sa1, label="COL_012")
#         col_xyz = dendropy.ColumnType(state_alphabet=sa2, label="COL_XYZ")
#         col_jkl = dendropy.ColumnType(state_alphabet=sa3, label="COL_JKL")
#         ca1.column_types = [col_012, col_xyz, col_jkl]
#
#         for t in taxon_set:
#             ca1[t] = dendropy.CharacterDataVector(self.get_standard_cells(col_012, "001122-??")) \
#                    + dendropy.CharacterDataVector(self.get_standard_cells(col_xyz, "XYZXYZ??-")) \
#                    + dendropy.CharacterDataVector(self.get_standard_cells(col_jkl, "JKJLKL-??")) \
#
#         for t, v in ca1.items():
#             self.assertEqual(len(v), 27)
#             for i, c in enumerate(v):
#                 if i >= 0 and i <= 8:
#                     self.assertTrue(c.column_type is col_012, [c.column_type, col_012])
#                     self.assertTrue(c.value in col_012.state_alphabet)
#                 elif i >= 9 and i <= 17:
#                     self.assertTrue(c.column_type is col_xyz, [c.column_type, col_xyz])
#                     self.assertTrue(c.value in col_xyz.state_alphabet)
#                 elif i >= 18 and i <= 26:
#                     self.assertTrue(c.column_type is col_jkl, [c.column_type, col_jkl])
#                     self.assertTrue(c.value in col_jkl.state_alphabet)
#
#         ca2 = dendropy.StandardCharacterArray(ca1)
#         self.assertTrue(ca1.taxon_set is ca2.taxon_set)
#         self.assertTrue(ca1.state_alphabets is not ca2.state_alphabets)
#         self.assertEqual(len(ca1.state_alphabets), len(ca2.state_alphabets))
#         self.assertEqual(len(ca1.state_alphabets), 3)
#         self.assertEqual(len(ca2.column_types), 3)
#         for ci, col1 in enumerate(ca1.column_types):
#             col2 = ca2.column_types[ci]
#             self.assertTrue(col2 is not col1)
#             self.assertTrue(col2.state_alphabet is not col1.state_alphabet)
#         for t, v1 in ca1.items():
#             v2 = ca2[t]
#             self.assertTrue(v1 is not v2)
#             for i, c1 in enumerate(v1):
#                 c2 = v2[i]
#                 self.assertTrue(c1 is not c2)
#                 self.assertTrue(c1.column_type is not c2.column_type)
#                 self.assertTrue(c1.column_type.state_alphabet is not c2.column_type.state_alphabet)
#
#     def get_standard_cells(self, col_type, symbols):
#         cells = [dendropy.CharacterDataCell(value=s, column_type=col_type) for s in col_type.state_alphabet.get_states(symbols=symbols)]
#         return cells
#
#     def get_standard_state_alphabet(self, symbols):
#         sa = dendropy.StateAlphabet()
#         for symbol in symbols:
#             sa.append(dendropy.StateAlphabetElement(symbol=symbol))
#         sa.append(dendropy.StateAlphabetElement(symbol="?",
#                                            multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
#                                            member_states=sa.get_states(symbols=symbols)))
#         sa.append(dendropy.StateAlphabetElement(symbol="-",
#                                            multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
#                                            member_states=sa.get_states(symbols=symbols)))
#         return sa
#
#
# class TreeInstantiationTest(unittest.TestCase):
#
#     def test_tree_init_from_newick(self):
#         newick_str = "((A,B),(C,D));"
#
#         # from file, using keywords
#         t1 = dendropy.Tree(istream=StringIO(newick_str), format="newick", oid="t1")
#         self.assertTrue(t1.oid == "t1", "'%s'" % t1.oid)
#         t1.debug_check_tree(_LOG)
#
#         # test copying
#         t2 = dendropy.Tree(t1)
#         t2.debug_check_tree(_LOG)
#         self.compare_tree_copies(t1, t2)
#
#         # from file, args
#         t3 = dendropy.Tree(StringIO(newick_str), "newick", taxon_set=t1.taxon_set)
#         t3.debug_check_tree(_LOG)
#         self.assertTrue(t3.taxon_set is t1.taxon_set)
#
#         # from file, mixed
#         t4 = dendropy.Tree(StringIO(newick_str), format="newick", taxon_set=t1.taxon_set)
#         t4.debug_check_tree(_LOG)
#         self.assertTrue(t4.taxon_set is t1.taxon_set)
#
#         # read from string
#         t5 = dendropy.Tree()
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
#         t1 = dendropy.Tree(istream=StringIO(nexus_str), format="nexus")
#         t2 = dendropy.Tree(StringIO(nexus_str), format="nexus")
#         t3 = dendropy.Tree(StringIO(nexus_str), "nexus")
#         t4 = dendropy.Tree()
#         t4.read_from_string(nexus_str, "nexus")
#         t5 = dendropy.Tree()
#         t5.read_from_file(StringIO(nexus_str), "nexus")
#         t6 = dendropy.Tree(StringIO(nexus_str), "nexus", taxon_set=t1.taxon_set)
#         self.assertTrue(t6.taxon_set is t1.taxon_set)
#         t7 = dendropy.Tree()
#         t7.read_from_string(nexus_str, "nexus", taxon_set=t1.taxon_set)
#         self.assertTrue(t7.taxon_set is t1.taxon_set)
#         for tx in (t1, t2, t3, t4, t5):
#             tx.debug_check_tree(_LOG)
#             tx.taxon_set.is_mutable = False
#             self.assertEqual(len(tx.taxon_set), 4, str([t.label for t in tx.taxon_set]))
#             self.assertTrue(tx.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
#         t8 = dendropy.Tree(t1)
#         self.compare_tree_copies(t1, t8)
#
#     def compare_tree_copies(self, t1, t2):
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
#
#     def test_treelist_init_from_newick(self):
#
#         newick_str = "((A,B),(C,D)); ((A,C),(B,D)); (A,(B,(C,D))); (A,(C,(B,D)));"
#
#         # from file, using keywords
#         tl1 = dendropy.TreeList(istream=StringIO(newick_str), format="newick", oid="t1")
#         self.assertTrue(tl1.oid == "t1", "'%s'" % tl1.oid)
#         self.assertEqual(len(tl1), 4)
#         self.assertEqual(len(tl1.taxon_set), 4)
#         for t in tl1:
#             t.debug_check_tree(_LOG)
#             self.assertTrue(tl1.taxon_set is t.taxon_set)
#
#         # test copying
#         tl2 = dendropy.TreeList(tl1)
#         self.assertTrue(tl2.taxon_set is tl1.taxon_set)
#         self.assertTrue(tl2.oid != tl1.oid)
#         self.assertTrue(tl2.label == tl1.label)
#         self.assertEqual(len(tl2.taxon_set), 4)
#         self.assertTrue(tl2.taxon_set.has_taxa(labels=["A", "B", "C", "D"]))
#         for ti, t1 in enumerate(tl1):
#             t2 = tl2[ti]
#             self.assertTrue(t2 is t1)
#
#         # from file, args
#         tl3 = dendropy.TreeList(StringIO(newick_str), "newick", taxon_set=tl1.taxon_set)
#         self.assertTrue(tl3.taxon_set is tl1.taxon_set)
#
#         # from file, mixed
#         tl4 = dendropy.TreeList(StringIO(newick_str), format="newick", taxon_set=tl1.taxon_set)
#         self.assertTrue(tl4.taxon_set is tl1.taxon_set)
#
#         # read from string
#         tl5 = dendropy.TreeList()
#         tl5.read_from_string(newick_str, format="newick")

if __name__ == "__main__":
    unittest.main()

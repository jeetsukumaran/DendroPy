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
NEWICK data read/write parse/format tests.
"""

import sys
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy import DataFormatError
import dendropy
from dendropy.dataio import newick

class NewickBasicParseTest(datatest.DataObjectVerificationTestCase):

    def testBadInit(self):
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c))a"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c)) (b,(a,c))"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c)) (d,(e,f))"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c)),"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c)))"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c)):"), format="NEWICK")
        self.assertRaises(DataFormatError, dendropy.TreeList, stream=StringIO("(a,(b,c))("), format="NEWICK")

    def testTreeListReaderDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        newick_str = datagen.reference_tree_list_newick_string()
        t_tree_list = dendropy.TreeList.get_from_string(newick_str, 'newick')
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None,
                ignore_taxon_order=True)

    def testTreeListReaderSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        newick_str = datagen.reference_tree_list_newick_string()
        t_tree_list = dendropy.TreeList.get_from_string(newick_str, 'newick', taxon_set=ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)

    def testReferenceTreeFileDistinctTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_path(pathmap.tree_source_path(datagen.reference_trees_filename(format="newick")), 'newick')
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None,
                ignore_taxon_order=True)

    def testReferenceTreeFileSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_path(pathmap.tree_source_path(datagen.reference_trees_filename(format="newick")),
                'newick',
                taxon_set=ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)
    def xtestToleratesExtraSemicolon(self):
        trees = dendropy.TreeList.get_from_string(
                """(((T1:1.1, T2:2.2)i1:4.0,(T3:3.3, T4:4.4)i2:4.0)i3:4.0,(T5:6.7, T6:7.2)i4:4.0)root:7.0;;""",
                "newick"
                )
        self.assertEquals(len(trees), 1)

class NewickEdgeLengthParsing(datatest.DataObjectVerificationTestCase):


    def testEdgeLengths1(self):
        trees = dendropy.TreeList.get_from_string(
                """(((T1:1.1, T2:2.2)i1:4.0,(T3:3.3, T4:4.4)i2:4.0)i3:4.0,(T5:6.7, T6:7.2)i4:4.0)root:7.0;""",
                "newick"
                )
        self.assertEquals(len(trees), 1)
        trees[0].debug_check_tree(self.logger)
        expected = {
            'T1': 1.1,
            'T2': 2.2,
            'i1': 4.0,
            'T3': 3.3,
            'T4': 4.4,
            'i2': 4.0,
            'i3': 4.0,
            'T5': 6.7,
            'T6': 7.2,
            'i4': 4.0,
            'root': 7.0
        }
        for nd in trees[0].postorder_node_iter():
            label = nd.taxon.label if nd.taxon is not None else nd.label
            self.assertAlmostEquals(nd.edge.length, expected[label])

    def testEdgeLengths2(self):
        trees = dendropy.TreeList.get_from_string("""
(((T1:1.242e-10, T2:213.31e-4)i1:3.44e-3,(T3:3.3e7, T4:4.4e+8)i2:4.0e+1)i3:4.0E-4,
(T5:6.7E+2, T6:7.2E-9)i4:4.0E8)root:7.0;
""", "newick")
        self.assertEquals(len(trees), 1)
        trees[0].debug_check_tree(self.logger)
        expected = {
            'T1': 1.242e-10,
            'T2': 213.31e-4,
            'i1': 3.44e-3,
            'T3': 3.3e7,
            'T4': 4.4e+8,
            'i2': 4.0e+1,
            'i3': 4.0e-4,
            'T5': 6.7e+2,
            'T6': 7.2e-9,
            'i4': 4.0e8,
            'root': 7.0
        }
        for nd in trees[0].postorder_node_iter():
            label = nd.taxon.label if nd.taxon is not None else nd.label
            self.assertAlmostEquals(nd.edge.length, expected[label])

    def testQuotedLabels(self):
        trees = dendropy.TreeList.get_from_string("""
((('T1 = 1.242e-10':1.242e-10,
'T2 is 213.31e-4':213.31e-4)i1:3.44e-3,
('T3 is a (nice) taxon':3.3e7,
T4:4.4e+8)'this is an internal node called "i2"':4.0e+1)i3:4.0E-4,
(T5:6.7E+2,
'and this so-called ''node'' is ("T6" with a length of ''7.2E-9'')':7.2E-9)i4:4.0E8)'this is the ''root\'\'\':7.0;
""", "newick")
        self.assertEquals(len(trees), 1)
        trees[0].debug_check_tree(self.logger)
        expected = {
            'T1 = 1.242e-10': 1.242e-10,
            'T2 is 213.31e-4': 213.31e-4,
            'i1': 3.44e-3,
            'T3 is a (nice) taxon': 3.3e7,
            'T4': 4.4e+8,
            'this is an internal node called "i2"': 4.0e+1,
            'i3': 4.0e-4,
            'T5': 6.7e+2,
            "and this so-called 'node' is (\"T6\" with a length of '7.2E-9')": 7.2e-9,
            'i4': 4.0e8,
            "this is the 'root'": 7.0
        }
        for nd in trees[0].postorder_node_iter():
            label = nd.taxon.label if nd.taxon is not None else nd.label
            self.assertAlmostEquals(nd.edge.length, expected[label])

class NewickTreeListWriterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()

    def testWriteTreeListDistinctTaxa(self):
        output_path = pathmap.named_output_path(filename="reference.trees.out.newick", suffix_timestamp=True)
        self.ref_tree_list.write_to_path(output_path, format="newick")
        t_tree_list = dendropy.TreeList.get_from_path(output_path, "newick")
        self.assertDistinctButEqualTreeList(
                self.ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None,
                ignore_taxon_order=True)

    def testWriteTreeListSameTaxa(self):
        output_path = pathmap.named_output_path(filename="reference.trees.out.newick", suffix_timestamp=True)
        self.ref_tree_list.write_to_path(output_path, format="newick")
        t_tree_list = dendropy.TreeList.get_from_path(output_path, "newick", taxon_set=self.ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                self.ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)

class NewickDocumentReaderTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        reference_tree_list = datagen.reference_tree_list()
        self.reference_dataset = dendropy.DataSet(reference_tree_list)

    def testBasicDocumentParseFromReader(self):
        reader = newick.NewickReader()
        test_dataset = reader.read(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(format="newick")))
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

    def testBasicDocumentParseFromRead(self):
        test_dataset = dendropy.DataSet()
        test_dataset.read(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(format="newick")), format="newick")
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

    def testBasicDocumentFromInit(self):
        test_dataset = dendropy.DataSet(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(format="newick")), format="newick")
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

class NewickDocumentWriterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        reference_tree_list = datagen.reference_tree_list()
        self.reference_dataset = dendropy.DataSet(reference_tree_list)

    def testRoundTrip(self):
        self.roundTripDataSetTest(self.reference_dataset, "newick", ignore_taxon_order=True)

if __name__ == "__main__":
    unittest.main()

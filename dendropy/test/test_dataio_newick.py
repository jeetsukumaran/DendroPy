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
NEWICK data read/write parse/format tests.
"""

import sys
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.utility.error import DataParseError
import dendropy
from dendropy.dataio import newick

class NewickBasicParseTest(datatest.DataObjectVerificationTestCase):

    def testBadInit(self):
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c))a"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c)) (b,(a,c))"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c)) (d,(e,f))"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c)),"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c)))"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c)):"), schema="NEWICK")
        self.assertRaises(DataParseError, dendropy.TreeList, stream=StringIO("(a,(b,c))("), schema="NEWICK")

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
        t_tree_list = dendropy.TreeList.get_from_path(pathmap.tree_source_path(datagen.reference_trees_filename(schema="newick")), 'newick')
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None,
                ignore_taxon_order=True)

    def testReferenceTreeFileSameTaxa(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList.get_from_path(pathmap.tree_source_path(datagen.reference_trees_filename(schema="newick")),
                'newick',
                taxon_set=ref_tree_list.taxon_set)
        self.assertDistinctButEqualTreeList(
                ref_tree_list,
                t_tree_list,
                distinct_taxa=False,
                equal_oids=None)
    def xtestToleratesExtraSemicolon(self):
        """Should an extra semicolon result in another (empty) tree being created ?
        MTH does not think so - but such strings are probably not legal newick,
        so whatever you want to do is OK with me
        """
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
            if nd.taxon is not None:
                label = nd.taxon.label
            else:
                label = nd.label
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
            if nd.taxon is not None:
                label = nd.taxon.label
            else:
                label = nd.label
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
            if nd.taxon is not None:
                label = nd.taxon.label
            else:
                label = nd.label
            self.assertAlmostEquals(nd.edge.length, expected[label])

class NewickTreeListWriterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.ref_tree_list = datagen.reference_tree_list()

    def testWriteTreeListDistinctTaxa(self):
        output_path = pathmap.named_output_path(filename="reference.trees.out.newick", suffix_timestamp=True)
        self.ref_tree_list.write_to_path(output_path, schema="newick")
        t_tree_list = dendropy.TreeList.get_from_path(output_path, "newick")
        self.assertDistinctButEqualTreeList(
                self.ref_tree_list,
                t_tree_list,
                distinct_taxa=True,
                equal_oids=None,
                ignore_taxon_order=True)

    def testWriteTreeListSameTaxa(self):
        output_path = pathmap.named_output_path(filename="reference.trees.out.newick", suffix_timestamp=True)
        self.ref_tree_list.write_to_path(output_path, schema="newick")
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
        test_dataset = reader.read(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(schema="newick")))
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

    def testBasicDocumentParseFromRead(self):
        test_dataset = dendropy.DataSet()
        test_dataset.read(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(schema="newick")), schema="newick")
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

    def testBasicDocumentFromInit(self):
        test_dataset = dendropy.DataSet(stream=pathmap.tree_source_stream(datagen.reference_trees_filename(schema="newick")), schema="newick")
        self.assertDistinctButEqual(self.reference_dataset, test_dataset, ignore_taxon_order=True)

class NewickDocumentWriterTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        reference_tree_list = datagen.reference_tree_list()
        self.reference_dataset = dendropy.DataSet(reference_tree_list)

    def testRoundTrip(self):
        self.roundTripDataSetTest(self.reference_dataset, "newick", ignore_taxon_order=True)

if __name__ == "__main__":
    unittest.main()

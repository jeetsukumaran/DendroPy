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

from dendropy.utility import messaging
from dendropy import tests
from dendropy.tests.test_tree_struct import TreeTraversalChecker
from dendropy.tests.rwtest import check_canonical_Pythonidae_cytb_tree_parse
import dendropy
from dendropy.dataio import newick

_LOG = messaging.get_logger(__name__)

def get_tree1():
    t1 = newick.read_tree_list(str="""
    (   (   A,
            (   C,
                E
            )D
        )B,
        (   (   H
            )I
        )G
    )F
    """)[0]
    return t1

def get_tree2():
    t2 = newick.read_tree_list(str="""
    (   (   (   (   T1,
                    (   (   T2,
                            T3
                        )i5,
                        T4
                    )i4
                )i3,
                (   T5,
                    T6,
                    T13
                )i6
            )i2,
            (   (   T7,
                    (   T8,
                        T9
                    )i9
                )i8,
                T10
            )i7
        )i1,
        T14,
        (   T11,
            T12
        )i10
    )i0;
    """)[0]
    return t2

class VerifyNewickParsedTreeTraversal(TreeTraversalChecker):

    def setUP(self):
        t1 = get_tree1()
        t2 = get_tree(2)
        t1.debug_check_tree(_LOG)
        t2.debug_check_tree(_LOG)
        tt = TreeTraversalChecker.setUp(self, t1, t2,
                "Tree Parsed From NEWICK Stream",
                set_oids_to_labels=True)

    def testTraversal(self):
        self.runChecks()

class VerifyNewickParsedTree(unittest.TestCase):

    def setUp(self):
        self.tree1 = get_tree1()
        self.tree2 = get_tree2()

    def testBasic(self):
        expected_tax_labels = ['A', 'C', 'H', 'E']
        self.assertEquals(len(self.tree1.taxon_set), 4)
        for t in self.tree1.taxon_set:
            self.assert_(t.label in expected_tax_labels)
        result_tax_labels = []
        for leaf in self.tree1.leaf_iter():
            self.assert_(leaf.taxon is not None)
            assert leaf.taxon.label in expected_tax_labels, \
                    "Unexpected taxon on tree: '%s'." % leaf.taxon.label
            result_tax_labels.append(leaf.taxon.label)
        for tax_label in expected_tax_labels:
            assert tax_label in result_tax_labels, \
                    "Expected taxon not found on tree: '%s'." % leaf.taxon.label

        expected_tax_labels = ["T%d" % i for i in xrange(1, 15)]
        self.assertEquals(len(self.tree2.taxon_set), 14)
        for t in self.tree2.taxon_set:
            assert t.label in expected_tax_labels
        result_tax_labels = []
        for leaf in self.tree2.leaf_iter():
            assert leaf.taxon is not None
            assert leaf.taxon.label in expected_tax_labels, \
                    "Unexpected taxon on tree: '%s'." % leaf.taxon.label
            result_tax_labels.append(leaf.taxon.label)
        for tax_label in expected_tax_labels:
            assert tax_label in result_tax_labels, \
                    "Expected taxon not found on tree: '%s'." % leaf.taxon.label

    def verifyTreesForGetTreesAndTreeIterTest(self,
            result,
            expected_tree_list=None,
            expected_taxa=None):
        self.assertEquals(len(result), 4)
        self.assertEquals(len(result.taxon_set), 3)
        for taxon in result.taxon_set:
            assert taxon.label in ['A', 'B', 'C']
        for tree in result:
            tree.debug_check_tree(_LOG)
            assert tree.taxon_set is result.taxon_set
            for nd in tree.leaf_iter():
                assert nd.taxon is not None
        if expected_tree_list is not None:
            assert result is expected_tree_list
        if expected_taxa is not None:
            assert result.taxon_set is expected_taxa

    def treestringForGetTreesAndTreeIterTest(self):
        return "(A,(B,C)); ((A,B),C); ((A,C),B); (A,B,C);"

    def testNewickGetTrees1(self):
        tc = newick.read_tree_list(str=self.treestringForGetTreesAndTreeIterTest())
        self.verifyTreesForGetTreesAndTreeIterTest(tc,
            expected_tree_list=None,
            expected_taxa=None)

    def testNewickGetTrees2(self):
        taxon_set = dendropy.TaxonSet()
        tc = newick.read_tree_list(str=self.treestringForGetTreesAndTreeIterTest(), taxon_set=taxon_set)
        self.verifyTreesForGetTreesAndTreeIterTest(tc,
            expected_tree_list=None,
            expected_taxa=taxon_set)

    def testNewickGetTrees2(self):
        tree_coll = dendropy.TreeList()
        tc = newick.read_tree_list(str=self.treestringForGetTreesAndTreeIterTest(), tree_list=tree_coll)
        self.verifyTreesForGetTreesAndTreeIterTest(tc,
            expected_tree_list=tree_coll,
            expected_taxa=tree_coll.taxon_set)

    def testNewickTreeIter1(self):
        tc = dendropy.TreeList()
        for t in newick.tree_source_iter(str=self.treestringForGetTreesAndTreeIterTest(), taxon_set=tc.taxon_set):
            tc.append(t)
        self.verifyTreesForGetTreesAndTreeIterTest(tc,
            expected_tree_list=None,
            expected_taxa=tc.taxon_set)

    def testNewickTreeIter2(self):
        tc = dendropy.TreeList()
        for t in newick.tree_source_iter(str=self.treestringForGetTreesAndTreeIterTest()):
            tc.append(t)
        self.verifyTreesForGetTreesAndTreeIterTest(tc,
            expected_tree_list=None,
            expected_taxa=tc.taxon_set)

    def testEdgeLengths1(self):
        trees = newick.read_tree_list(str="""
(((T1:1.1, T2:2.2)i1:4.0,(T3:3.3, T4:4.4)i2:4.0)i3:4.0,(T5:6.7, T6:7.2)i4:4.0)root:7.0;
""")
        assert len(trees) == 1
        trees[0].debug_check_tree(_LOG)
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
        trees = newick.read_tree_list(str="""
(((T1:1.242e-10, T2:213.31e-4)i1:3.44e-3,(T3:3.3e7, T4:4.4e+8)i2:4.0e+1)i3:4.0E-4,
(T5:6.7E+2, T6:7.2E-9)i4:4.0E8)root:7.0;
""")
        assert len(trees) == 1
        trees[0].debug_check_tree(_LOG)
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
        trees = newick.read_tree_list(str="""
((('T1 = 1.242e-10':1.242e-10,
'T2 is 213.31e-4':213.31e-4)i1:3.44e-3,
('T3 is a (nice) taxon':3.3e7,
T4:4.4e+8)'this is an internal node called "i2"':4.0e+1)i3:4.0E-4,
(T5:6.7E+2,
'and this so-called ''node'' is ("T6" with a length of ''7.2E-9'')':7.2E-9)i4:4.0E8)'this is the ''root''':7.0;
""")
        assert len(trees) == 1
        trees[0].debug_check_tree(_LOG)
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

class NewickWriterTest(unittest.TestCase):

    def testReadTreeList(self):
        check_canonical_Pythonidae_cytb_tree_parse(
                reader = newick.NewickReader(),
                srcpath=tests.data_source_path("pythonidae_cytb.newick.tre"),
                logger=_LOG,
                underscore_substitution=True)

    def testWriteTreeList(self):
        _LOG.info("Reading in trees for NEWICK writing test")
        reader = newick.NewickReader()
        ds1 = reader.read(file=open(tests.data_source_path("pythonidae_cytb.newick.tre"), "rU"))

        outfile = tempfile.NamedTemporaryFile()
        _LOG.info("Writing trees to temporary file '%s'" % outfile.name)
        writer = newick.NewickWriter(dataset=ds1)
        writer.write(file=outfile)
        outfile.flush()

        _LOG.info("Re-reading trees")
        check_canonical_Pythonidae_cytb_tree_parse(
                reader = newick.NewickReader(),
                srcpath=outfile.name,
                logger=_LOG,
                underscore_substitution=True)

if __name__ == "__main__":
    unittest.main()

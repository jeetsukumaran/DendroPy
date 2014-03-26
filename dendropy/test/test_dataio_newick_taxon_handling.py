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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import dendropy

class NewickTreeTaxonHandlingTest(unittest.TestCase):

    def test_basic_taxa(self):
        s = "(a1:3.14e-2,(b2:1.2,(c3:0.5,d4:0.7)e5:111)f6:222)g7:333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a1": 3.14e-2,
                "b2": 1.2,
                "c3": 0.5,
                "d4": 0.7,
                "e5": 111,
                "f6": 222,
                "g7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

    def test_quoted_underscores(self):
        s = "('a_1':3.14e-2,('b_2':1.2,('c_3':0.5,'d_4':0.7)'e_5':111)'f_6':222)'g_7':333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a_1": 3.14e-2,
                "b_2": 1.2,
                "c_3": 0.5,
                "d_4": 0.7,
                "e_5": 111,
                "f_6": 222,
                "g_7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

    def test_unquoted_underscores(self):
        s = "(a_1:3.14e-2,(b_2:1.2,(c_3:0.5,d_4:0.7)e_5:111)f_6:222)g_7:333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a 1": 3.14e-2,
                "b 2": 1.2,
                "c 3": 0.5,
                "d 4": 0.7,
                "e 5": 111,
                "f 6": 222,
                "g 7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])


if __name__ == "__main__":
    unittest.main()

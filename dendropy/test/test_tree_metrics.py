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
Tests of tree metrics.
"""

import random
import unittest
import math
from cStringIO import StringIO

import dendropy
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import extendedtest

class TreePHGammTest(unittest.TestCase):

    def runTest(self):
        newick_str = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        tree = dendropy.Tree(stream=StringIO(newick_str), schema="newick")
        self.assertAlmostEqual(tree.pybus_harvey_gamma(tree), 0.546276, 4)

class TreeCompareTests(extendedtest.ExtendedTestCase):

    def setUp(self):
        self.tree_list1 = datagen.reference_tree_list()
        self.tree_list2 = datagen.reference_tree_list()

    def testNonMutatingDistinctTaxonSetSameStructComparisons(self):
        tl1_ts = self.tree_list1.taxon_set
        tl2_ts = self.tree_list2.taxon_set
        self.assertIsNot(tl1_ts, tl2_ts)
        for i, t1 in enumerate(self.tree_list1):
            t2 = self.tree_list2[i]
            t1_ts = t1.taxon_set
            t2_ts = t2.taxon_set
            self.assertIsNot(t1_ts, t2_ts)
            self.assertEqual(t1.symmetric_difference(t2), 0)
            self.assertAlmostEqual(t1.euclidean_distance(t2), 0)
            self.assertAlmostEqual(t1.robinson_foulds_distance(t2), 0)
            self.assertIs(t1.taxon_set, t1_ts)
            self.assertIs(t2.taxon_set, t2_ts)
        self.assertIs(self.tree_list1.taxon_set, tl1_ts)
        self.assertIs(self.tree_list2.taxon_set, tl2_ts)

class FrequencyOfSplitsTest(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(
                src=pathmap.tree_source_path('pythonidae.random.bd0301.tre'),
                schema='nexus')

    def testCount1(self):
        split_leaves = ['Python regius', 'Apodora papuana']
        f = self.trees.frequency_of_split(labels=split_leaves)
        self.assertAlmostEqual(f, 0.04)

    def testRaisesIndexError(self):
        split_leaves = ['Bad Taxon', 'Apodora papuana']
        self.assertRaises(IndexError, self.trees.frequency_of_split, labels=split_leaves)

if __name__ == "__main__":
    unittest.main()

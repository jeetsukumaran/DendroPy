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

class TreePHGammTest(unittest.TestCase):

    def runTest(self):
        newick_str = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        tree = dendropy.Tree(stream=StringIO(newick_str), schema="newick")
        self.assertAlmostEqual(tree.pybus_harvey_gamma(tree), 0.546276, 4)

class FrequencyOfSplitsTest(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(
                src=pathmap.tree_source_path('pythonidae.random.tre'),
                schema='nexus')

    def testCount1(self):
        split_leaves = ['Python regius', 'Apodora papuana']
        f = self.trees.frequency_of_split(labels=split_leaves)
        self.assertAlmostEqual(f, 0.02)

    def testRaisesIndexError(self):
        split_leaves = ['Bad Taxon', 'Apodora papuana']
        self.assertRaises(IndexError, self.trees.frequency_of_split, labels=split_leaves)

if __name__ == "__main__":
    unittest.main()

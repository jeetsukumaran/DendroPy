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
Tests the interaction with APE.
"""

import unittest
from dendropy.utility import error
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.test.support import pathmap
from dendropy.interop import ape
import dendropy

class DataRoundTrip(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.trees = datagen.reference_tree_list()
        self.dna_chars = datagen.reference_dna_array()
        self.std_chars = datagen.reference_standard_array()

    def testTreeRoundTrip(self):
        for i, t1 in enumerate(self.trees):
            ape_t = ape.as_ape_object(t1)
            rt = ape.as_dendropy_object(ape_t, taxon_set=t1.taxon_set)
            self.assertTrue(isinstance(rt, dendropy.Tree), type(rt))
            self.assertDistinctButEqual(t1, rt, distinct_taxa=False)

    def testTreeListRoundTrip(self):
        ape_t = ape.as_ape_object(self.trees)
        rt = ape.as_dendropy_object(ape_t, taxon_set=self.trees.taxon_set)
        self.assertTrue(isinstance(rt, dendropy.TreeList), type(rt))
        self.assertDistinctButEqual(self.trees, rt, distinct_taxa=False)

    def testDnaRoundTrip(self):
        ape_c = ape.as_ape_object(self.dna_chars)
        dp_c = ape.as_dendropy_object(ape_c, taxon_set=self.dna_chars.taxon_set)
        self.assertDistinctButEqual(self.dna_chars, dp_c, distinct_taxa=False)

#    def testStandardRoundTrip(self):
#        ape_c = ape.as_ape_object(self.std_chars)
#        st_c = ape.as_dendropy_object(ape_c, taxon_set=self.std_chars.taxon_set)
#        self.assertDistinctButEqual(self.std_chars, st_c, distinct_taxa=False)

class ApeFunctions(unittest.TestCase):

    def setUp(self):
        self.tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('single-ultrametric.tre'), 'nexus')

    def testBirthDeath(self):
        bd = ape.birthdeath(self.tree)
        expected = {
            'deviance': -168.44349668256115,
            'log-likelihood': -84.221748341280573,
            'd/b': 0.0,
            'd/b s.e.': 0.0,
            'b-d': 934.02440046546621,
            'b-d s.e.': 199.15458447446244,
        }
        for k, v in expected.items():
            self.assertAlmostEqual(v, bd[k])

    def testBirthDeathExt(self):
        t = dendropy.Tree(self.tree)
        num_species = [9, 7, 4, 5, 9, 6, 7, 2, 8, 6, 3, 9, 6]
        for i, n in enumerate(t.leaf_iter()):
            n.num_species = num_species[i]
        bd = ape.bd_ext(t)
#        expected = {
#            'deviance': -44.609589999999997,
#            'log-likelihood': 22.3048,
#            'd/b': 0.99938660000000001,
#            'd/b s.e.': 0.00011987800000000001,
#            'b-d': 1.560403,
#            'b-d s.e.': 0.38816460000000003,
#        }
        expected = {
            'deviance': -140.75389999999999,
            'log-likelihood': 70.376930000000002,
            'd/b': 0.0,
            'd/b s.e.': 0.0,
            'b-d': 2900.3229999999999,
            'b-d s.e.': 273.3707,
        }
        for k, v in expected.items():
            self.assertAlmostEqual(v, bd[k])

if __name__ == "__main__":
    unittest.main()


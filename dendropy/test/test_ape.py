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
Tests the interaction with APE.
"""

import unittest
from dendropy.utility import error
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.test.support import pathmap
from dendropy.interop import ape
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

import dendropy

if not ape.DENDROPY_APE_INTEROPERABILITY:
    _LOG.warn("APE interoperability not available: skipping APE tests")
else:
    class DataRoundTrip(datatest.DataObjectVerificationTestCase):

        def setUp(self):
            self.trees = datagen.reference_tree_list()
            self.dna_chars = datagen.reference_dna_matrix()
            self.std_chars = datagen.reference_standard_matrix()

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
            self.tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('bird_orders.nex'), 'nexus')

        def testBirthDeath(self):
            bd = ape.birthdeath(self.tree)
            expected = {
                'deviance': 81.20744,
                'log-likelihood': -40.60372,
                'd/b': 0.0,
                'd/b s.e.': 0.0,
                'b-d': 0.03909895,
                'b-d s.e.': 0.006048519,
            }
            for k, v in expected.items():
                self.assertAlmostEqual(v, bd[k], 2)

        def testBirthDeathExt(self):
            t = dendropy.Tree(self.tree)
            num_species = [10, 47, 69, 214, 161, 17,
                355, 51, 56, 10, 39, 152,
                6, 143, 358, 103, 319,
                23, 291, 313, 196, 1027, 5712]
            for i, n in enumerate(t.taxon_set):
                n.num_species = num_species[i]
            bd = ape.bd_ext(t)
            expected = {
                'deviance': 289.1639,
                'log-likelihood': -144.5820,
                'd/b': 0.0,
                'd/b s.e.': 0.0,
                'b-d': 0.2866789,
                'b-d s.e.': 0.007215593,
            }
            for k, v in expected.items():
                self.assertAlmostEqual(v, bd[k], 3)

    if __name__ == "__main__":
        unittest.main()


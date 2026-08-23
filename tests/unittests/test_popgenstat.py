#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Tests population genetic statistic calculation.
"""

import unittest
import math
import dendropy
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest
from support import pathmap
from dendropy.utility import messaging
from dendropy.calculate import popgenstat
_LOG = messaging.get_logger(__name__)

class TajimasDTests(dendropytest.ExtendedTestCase):

    def setUp(self):
        s = """\
            >s1
            ATAATAAAAA AATAATAAAA AAATAAAAAA AATAAAAAAA A
            >s2
            AAAAAAAATA AATAATAAAA AAATAAAAAA AAAAAAAAAA A
            >s3
            AAAATAAAAA TATAATAAAA AAATATAAAA AAAAAAAAAA A
            >s4
            AAAAAAAAAA AATAATAAAA AAATAAATAA ATAAAAAAAA A
            >s5
            AAAATAAAAA AAATATAAAA AAATAAAAAA AAAAAAAAAA A
            >s6
            AAAATAAAAA AAAAATAAAA AAAAAAAAAA AAAAATAAAA A
            >s7
            AAAAAATAAA AATAATAAAA AAATAAAAAA AAAAAAAAAA A
            >s8
            AAAAAAAAAA AAAAATAAAA AAATAAAAAA AAAAAAAAAT A
            >s9
            AAAAAAAAAA AAAAAAAAAA AAATAAAAAA AAAAAAAAAA A
            >s10
            AAAAAAAAAA AAAAATAAAA AAATAATAAA AAAAAAAAAA A"""
        self.matrix = dendropy.DnaCharacterMatrix.get_from_string(s, 'fasta')

        s_with_missing = """\
            >s1
            ATAATAAAAA AATAATAAAA AAATAAAAAA AATAAAAAAA A
            >s2
            AAAAAAAATA AATAATAAAA AAATAAAAAA AAAAAAAAAA A
            >s3
            AAAATAAAAA TATAATAAAA AAATATAAAA ?AAAAAAAAA ?
            >s4
            AAAAAAAAAA AATAATAAAA AAATAAATAA ATAAAAAAAA A
            >s5
            AAAATAAAAA AAATATAAAA AAATAAAAAA AAAAAAAAAA A
            >s6
            ?AAATAAAAA AAAAATAAAA AAAAAAAAAA AAAAATAAAA A
            >s7
            AAAAAATAAA AATAATAAAA AAATAAAAAA AAAAAAAAAA A
            >s8
            AAAAAAAAAA AAAAATAAAA AAATAAAAAA AAAAAAAAAT A
            >s9
            AAAAAAAAAA AAAAAAAAAA AAATAAAAAA ?AAAAAAAAA A
            >s10
            AAAAAAAAAA AAAAATAAAA AAATAATAAA AAAAAAAAAA A"""
        self.matrix_with_missing = dendropy.DnaCharacterMatrix.get_from_string(s_with_missing, 'fasta')

    def testTajimasD(self):
        self.assertAlmostEqual(popgenstat.tajimas_d(self.matrix), -1.44617198561, 4)

    def testTajimasD_with_missing(self):
        self.assertAlmostEqual(popgenstat.tajimas_d(self.matrix_with_missing), -1.44617198561, 4)

class SinglePopTest(dendropytest.ExtendedTestCase):

    data = dendropy.DnaCharacterMatrix.get_from_path(pathmap.char_source_path('COII_Apes.nex'), schema="nexus")

    def test_num_segregating_sites(self):
        self.assertEqual(popgenstat.num_segregating_sites(self.data, ignore_uncertain=True), 183)

    def test_average_number_of_pairwise_differences(self):
        self.assertAlmostEqual(popgenstat.average_number_of_pairwise_differences(self.data, ignore_uncertain=True),  62.75000, 4)

    def test_nucleotide_diversity(self):
        self.assertAlmostEqual(popgenstat.nucleotide_diversity(self.data, ignore_uncertain=True), 0.09174, 4)

    def test_tajimas_d(self):
        self.assertAlmostEqual(popgenstat.tajimas_d(self.data, ignore_uncertain=True), 1.12467, 4)

    def test_wattersons_theta(self):
        self.assertAlmostEqual(popgenstat.wattersons_theta(self.data, ignore_uncertain=True), 49.00528, 4)

class PopulationPairSummaryStatisticsTests(dendropytest.ExtendedTestCase):

    def testPopulationPairSummaryStatistics(self):
        seqs = dendropy.DnaCharacterMatrix.get_from_path(pathmap.char_source_path('orti.nex'), schema="nexus")
        p1 = []
        p2 = []
        for idx, t in enumerate(seqs.taxon_namespace):
            if t.label.startswith('EPAC'):
                p1.append(seqs[t])
            else:
                p2.append(seqs[t])
        pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)
        self.assertAlmostEqual(pp.average_number_of_pairwise_differences, 11.28063, 4)
        self.assertAlmostEqual(pp.average_number_of_pairwise_differences_between, 16.119047619, 4)
        self.assertAlmostEqual(pp.average_number_of_pairwise_differences_within, 10.2191697192, 4)
        self.assertAlmostEqual(pp.average_number_of_pairwise_differences_net, 5.89987789988, 4)
        self.assertEqual(pp.num_segregating_sites, 29)
        self.assertAlmostEqual(pp.wattersons_theta, 7.85734688643, 4)
        self.assertAlmostEqual(pp.tajimas_d, 1.65318627677, 4)
        self.assertAlmostEqual(pp.wakeleys_psi, 0.8034976, 2)

class PopulationPairAmbiguityTests(dendropytest.ExtendedTestCase):
    "Ambiguous / gap / no-data bases are treated as missing between populations."

    def _dxy(self, pop2_char):
        # Both populations carry a within-population difference (so calc() does
        # not divide by zero); only pop2 sequence-a site 3 varies between calls.
        s = ">P1_a\nCAAAA\n>P1_b\nAAAAA\n>P2_a\nAA%sAG\n>P2_b\nAAAAA\n" % pop2_char
        matrix = dendropy.DnaCharacterMatrix.get_from_string(s, "fasta")
        p1 = [matrix[t] for t in matrix.taxon_namespace if t.label.startswith("P1")]
        p2 = [matrix[t] for t in matrix.taxon_namespace if t.label.startswith("P2")]
        return popgenstat.PopulationPairSummaryStatistics(
                p1, p2).average_number_of_pairwise_differences_between

    def test_ambiguous_and_gap_treated_as_missing(self):
        # '?' excludes the varying site, leaving d_xy = 1.0; 'N' and '-' must match.
        self.assertAlmostEqual(self._dxy("?"), 1.0, 8)
        self.assertAlmostEqual(self._dxy("N"), 1.0, 8)
        self.assertAlmostEqual(self._dxy("-"), 1.0, 8)

    def test_partial_ambiguity_still_compared(self):
        # A partial code (R = A|G) is not missing, so it must still be compared.
        self.assertNotAlmostEqual(self._dxy("R"), self._dxy("?"), 8)

if __name__ == "__main__":
    unittest.main()

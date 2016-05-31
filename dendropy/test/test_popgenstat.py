#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests population genetic statistic calculation.
"""

import unittest
import math
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap
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

if __name__ == "__main__":
    unittest.main()

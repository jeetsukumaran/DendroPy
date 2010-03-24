#! /usr/bin/env python

############################################################################
##  test_tree_io.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
Tests population genetic statistic calculation.
"""

import unittest
import StringIO
import math
import dendropy
from dendropy.test.support import extendedtest
from dendropy.test.support import pathmap
from dendropy.utility import messaging
from dendropy import popgenstat
_LOG = messaging.get_logger(__name__)

class TajimasDTests(extendedtest.ExtendedTestCase):

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

    def testTajimasD(self):
        self.assertAlmostEqual(popgenstat.tajimas_d(self.matrix), -1.44617198561, 4)

class SinglePopTest(extendedtest.ExtendedTestCase):

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

class PopulationPairSummaryStatisticsTests(extendedtest.ExtendedTestCase):

    def testPopulationPairSummaryStatistics(self):
        seqs = dendropy.DnaCharacterMatrix.get_from_path(pathmap.char_source_path('orti.nex'), schema="nexus")
        p1 = []
        p2 = []
        for idx, t in enumerate(seqs.taxon_set):
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


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
Test the PAUP* wrapper
"""

import os
import sys
import csv

import unittest
from dendropy.test.support import pathmap
from dendropy.test.support.framework import ExtendedTestCase
from dendropy.utility import containers
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

from dendropy import splitcalc
from dendropy.utility import paup

if "PAUP_PATH" in os.environ:
    PAUP_PATH = os.environ["PAUP_PATH"]
else:
    PAUP_PATH = "paup"

class PaupWrapperRepToSplitMaskTest(unittest.TestCase):

    def testUnnormalized(self):
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, ".", "*")[::-1]
            r = paup.paup_group_to_mask(s, normalized=False)
            self.assertEqual(r, i, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(r, 8)))

    def testNormalized0(self):
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, "*", ".")[::-1]
            r = paup.paup_group_to_mask(s, normalized=True)
            normalized = containers.NormalizedBitmaskDict.normalize(i, 0xFF)
            self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(normalized, 8)))

    def testNormalized1(self):
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, ".", "*")[::-1]
            r = paup.paup_group_to_mask(s, normalized=True)
            normalized = containers.NormalizedBitmaskDict.normalize(i, 0xFF)
            self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(normalized, 8)))

class PaupWrapperSplitsParse(ExtendedTestCase):

    def setUp(self):
        self.tree_filepath = None
        self.taxa_filepath = None
        self.splitscsv_filepath = None
        self.expected_num_trees = None
        self.expected_split_freqs = None

    def populate_test(self, tree_filepath, splitscsv_filepath, expected_num_trees):
        self.tree_filepath = pathmap.data_source_path(tree_filepath)
        self.taxa_filepath = self.tree_filepath
        self.splitscsv_filepath = pathmap.data_source_path(splitscsv_filepath)
        self.expected_num_trees = expected_num_trees
        self.expected_split_freqs = dict([ (s[0], int(s[1])) for s in csv.reader(open(self.splitscsv_filepath, "rU"))])

    def count_splits(self, unrooted=True):
        if self.tree_filepath is None:
            _LOG.warning("Null Test Case")
            return
        p = paup.PaupRunner()
        p.stage_execute_file(self.taxa_filepath, clear_trees=True)
        p.stage_list_taxa()
        p.stage_load_trees(tree_filepaths=[self.tree_filepath], unrooted=unrooted)
        p.stage_count_splits()
        p.run()
        taxon_set = p.parse_taxon_set()
        tree_count, bipartition_counts = p.parse_group_freqs()

        self.assertEqual(self.expected_num_trees, tree_count)
        self.assertEqual(len(self.expected_split_freqs), len(bipartition_counts))
        for g in self.expected_split_freqs:
            self.assertIsContainedIn(g, bipartition_counts)
            self.assertEqual(self.expected_split_freqs[g], bipartition_counts[g])

        sd = paup.build_split_distribution(bipartition_counts,
                                           tree_count,
                                           taxon_set,
                                           unrooted=unrooted)
        sf = sd.split_frequencies
        for g in bipartition_counts:
            s = paup.paup_group_to_mask(g, normalized=unrooted)
            self.assertIsContainedIn(s, sd.splits)
            self.assertIsContainedIn(s, sd.split_counts)
            self.assertEqual(sd.split_counts[s], bipartition_counts[g])
            self.assertEqual(sd.total_trees_counted, self.expected_num_trees)
            self.assertAlmostEqual(sf[s], float(bipartition_counts[g]) / self.expected_num_trees)

class PaupWrapperSplitsParseTest1(PaupWrapperSplitsParse):

    def setUp(self):
        self.populate_test(["trees","feb032009.tre"],
                                ["trees", "feb032009.splits.csv"],
                                100)

    def runTest(self):
        self.count_splits()

class PaupWrapperTaxaParse(ExtendedTestCase):

    def setUp(self):
        self.taxa_filepath = None
        self.expected_taxlabels = None

    def check_labels(self):
        p = paup.PaupRunner()
        p.stage_execute_file(self.taxa_filepath)
        p.stage_list_taxa()
        p.run()
        taxon_set = p.parse_taxon_set()
        self.assertEqual(len(taxon_set), len(self.expected_taxlabels))
        for i, t in enumerate(taxon_set):
            self.assertEqual(t.label, self.expected_taxlabels[i])

class PaupWrapperTaxaParseTest1(PaupWrapperTaxaParse):

    def setUp(self):
        self.taxa_filepath = pathmap.data_source_path(["trees", "feb032009.tre"])
        self.expected_taxlabels = ("T01", "T02", "T03", "T04", "T05", "T06",
                                   "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14",
                                   "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22",
                                   "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
                                   "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38",
                                   "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46",
                                   "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54",
                                   "T55", "T56", "T57", "T58", "T59")

    def runTest(self):
        self.check_labels

class PaupWrapperTaxaParseTest2(PaupWrapperTaxaParse):

    def setUp(self):
        self.taxa_filepath = pathmap.data_source_path(["chars", "primates.chars.nexus"])
        self.expected_taxlabels = ("Lemur catta", "Homo sapiens",
                "Pan", "Gorilla", "Pongo", "Hylobates", "Macaca fuscata",
                "Macaca mulatta", "Macaca fascicularis", "Macaca sylvanus",
                "Saimiri sciureus", "Tarsius syrichta", )

    def runTest(self):
        self.check_labels

if __name__ == "__main__":
    unittest.main()

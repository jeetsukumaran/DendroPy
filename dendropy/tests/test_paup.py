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
import subprocess
import tempfile
import re
import csv

import unittest
import dendropy.tests
from dendropy.utility import containers
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

from dendropy import splitcalc
from dendropy.utility import paup

if "PAUP_PATH" in os.environ:
    PAUP_PATH = os.environ["PAUP_PATH"]
else:
    PAUP_PATH = "paup"

class PaupWrapperDumbTests(unittest.TestCase):
    """
    Checks basic running of PAUP commands, and correct parsing/extraction of
    output.
    """

    class SplitsTestCase(object):

        def __init__(self, tree_filepath, taxa_filepath, splitscsv_filepath, num_trees):
            self.tree_filepath = tree_filepath
            self.taxa_filepath = taxa_filepath
            self.splitscsv_filepath = splitscsv_filepath
            self.num_trees = num_trees

        def _get_treefilepath(self):
            return self._tree_filepath
        def _set_treefilepath(self, f):
            self._tree_filepath = dendropy.tests.data_source_path(f)
        tree_filepath = property(_get_treefilepath, _set_treefilepath)

        def _get_taxafilepath(self):
            return self._taxa_filepath
        def _set_taxafilepath(self, f):
            if f is not None:
                self._taxa_filepath = dendropy.tests.data_source_path(f)
            else:
                self._taxa_filepath = None
        taxa_filepath = property(_get_taxafilepath, _set_taxafilepath)

        def _get_splitscsvfilepath(self):
            return self._splitscsv_filepath
        def _set_splitscsvfilepath(self, f):
            self._splitscsv_filepath = dendropy.tests.data_source_path(f)
        splitscsv_filepath = property(_get_splitscsvfilepath, _set_splitscsvfilepath)

    def setUp(self):

        self.splits_test_cases = [
            PaupWrapperDumbTests.SplitsTestCase(
                ["trees","feb032009.tre"],
                ["trees", "feb032009.tre"],
                ["trees", "feb032009.splits.csv"],
                100),
            ]
        self.taxa_test_cases = (
            (["trees", "feb032009.tre"],
                ("T01", "T02", "T03", "T04", "T05", "T06",
                "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14",
                "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22",
                "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
                "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38",
                "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46",
                "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54",
                "T55", "T56", "T57", "T58", "T59")),
            (["chars", "primates.chars.nexus"],
                ("Lemur catta", "Homo sapiens",
                "Pan", "Gorilla", "Pongo", "Hylobates", "Macaca fuscata",
                "Macaca mulatta", "Macaca fascicularis", "Macaca sylvanus",
                "Saimiri sciureus", "Tarsius syrichta", ))
        )

    def check_taxon_set(self, filename, taxlabels):
        """Loads a taxa block from `filename`, make sure taxa returned match
        `taxlabels`"""
        p = paup.PaupRunner()
        p.stage_execute_file(dendropy.tests.data_source_path(filename))
        p.stage_list_taxa()
        p.run()
        taxon_set = p.parse_taxon_set()
        assert len(taxon_set) == len(taxlabels)
        for i, t in enumerate(taxon_set):
            assert t.label == taxlabels[i]

    def check_group_freqs(self, treefile, taxafile, exp_tree_count, group_freqs, unrooted=True):
        """Calculates group frequencies from `filename`, make sure that
        frequencies match `group_freqs` (given as dictionary of PAUP* group
        strings and their counts for the file)."""
        p = paup.PaupRunner()
        p.stage_execute_file(taxafile, clear_trees=True)
        p.stage_list_taxa()
        p.stage_load_trees(tree_filepaths=[treefile], unrooted=unrooted)
        p.stage_count_splits()
        p.run()
        taxon_set = p.parse_taxon_set()
        tree_count, bipartition_counts = p.parse_group_freqs()

        assert exp_tree_count == tree_count, "%s != %s" % (exp_tree_count, tree_count)
        assert len(group_freqs) == len(bipartition_counts), "%d != %d" % (len(group_freqs), len(bipartition_counts))
        for g in group_freqs:
            assert g in bipartition_counts
            assert group_freqs[g] == bipartition_counts[g], \
                "%s != %s" % (group_freqs[g], bipartition_counts[g])

        sd = paup.build_split_distribution(bipartition_counts,
                                           tree_count,
                                           taxon_set,
                                           unrooted=unrooted)

        sf = sd.split_frequencies
        for g in bipartition_counts:
            s = paup.paup_group_to_mask(g, normalized=unrooted)
            assert s in sd.splits
            assert s in sd.split_counts
            assert sd.split_counts[s] == bipartition_counts[g]
            assert sd.total_trees_counted == exp_tree_count
            self.assertAlmostEqual(sf[s], float(bipartition_counts[g]) / exp_tree_count)

    def testGroupRepToSplitMask(self):
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, ".", "*")[::-1]
            r = paup.paup_group_to_mask(s, normalized=False)
            assert r == i, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(r, 8))
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, "*", ".")[::-1]
            r = paup.paup_group_to_mask(s, normalized=True)
            normalized = containers.NormalizedBitmaskDict.normalize(i, 0xFF)
            assert r == normalized, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(normalized, 8))
        for i in xrange(0xFF):
            s = splitcalc.split_as_string(i, 8, ".", "*")[::-1]
            r = paup.paup_group_to_mask(s, normalized=True)
            normalized = containers.NormalizedBitmaskDict.normalize(i, 0xFF)
            assert r == normalized, "%s  =>  %s  =>  %s" \
                % (splitcalc.split_as_string(i, 8), s, splitcalc.split_as_string(normalized, 8))

    def testTaxaBlock(self):
        for i in self.taxa_test_cases:
            self.check_taxon_set(i[0], i[1])

    def testGroupFreqs(self):
        for stc in self.splits_test_cases:
            splits = dict([ (s[0], int(s[1])) for s in csv.reader(open(stc.splitscsv_filepath, "rU"))])
            _LOG.info("Checking splits in %s" % stc.tree_filepath)
            self.check_group_freqs(stc.tree_filepath,
                stc.taxa_filepath, stc.num_trees, splits)

if __name__ == "__main__":
    unittest.main()

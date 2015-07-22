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
Test the PAUP* wrapper
"""

import os
import sys
import csv
import collections
import unittest

from dendropy.test.support import pathmap
from dendropy.test.support import paupsplitsreference
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.utility import messaging
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
_LOG = messaging.get_logger(__name__)

from dendropy.utility import bitprocessing
from dendropy.interop import paup
from dendropy import Bipartition

if not paup.DENDROPY_PAUP_INTEROPERABILITY:
    _LOG.warn("PAUP interoperability not available: skipping PAUP tests")
else:

    class PaupWrapperRepToSplitMaskTest(unittest.TestCase):

        def setUp(self):
            self.ps = paup.PaupService()

        def testUnnormalized(self):
            for i in range(0xFF):
                s = bitprocessing.int_as_bitstring(i, 8, ".", "*")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=False)
                self.assertEqual(r, i, "%s  =>  %s  =>  %s" \
                    % (bitprocessing.int_as_bitstring(i, 8), s, bitprocessing.int_as_bitstring(r, 8)))

        def testNormalized0(self):
            for i in range(0xFF):
                s = bitprocessing.int_as_bitstring(i, 8, "*", ".")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=True)
                normalized = Bipartition.normalize_bitmask(i, 0xFF, 1)
                self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                    % (bitprocessing.int_as_bitstring(i, 8), s, bitprocessing.int_as_bitstring(normalized, 8)))

        def testNormalized1(self):
            for i in range(0xFF):
                s = bitprocessing.int_as_bitstring(i, 8, ".", "*")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=True)
                normalized = Bipartition.normalize_bitmask(i, 0xFF, 1)
                self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                    % (bitprocessing.int_as_bitstring(i, 8), s, bitprocessing.int_as_bitstring(normalized, 8)))

    class PaupWrapperSplitsParse(ExtendedTestCase):

        def check_splits_counting(self,
                tree_filename,
                taxa_definition_filepath,
                splits_filename,
                paup_as_rooted,
                paup_use_tree_weights,
                paup_burnin,
                expected_taxon_labels,
                expected_is_rooted,
                expected_num_trees,
                ):
            tree_filepath = pathmap.tree_source_path(tree_filename)
            paup_service = paup.PaupService()
            result = paup_service.count_splits_from_files(
                    tree_filepaths=[tree_filepath],
                    taxa_definition_filepath=taxa_definition_filepath,
                    is_rooted=paup_as_rooted,
                    use_tree_weights=paup_use_tree_weights,
                    burnin=paup_burnin,
                    )
            num_trees = result["num_trees"]
            bipartition_counts = result["bipartition_counts"]
            bipartition_freqs = result["bipartition_freqs"]
            taxon_namespace = result["taxon_namespace"]
            is_rooted = result["is_rooted"]

            # check taxon namespace
            self.assertEqual(len(taxon_namespace), len(expected_taxon_labels))
            for taxon, expected_label in zip(taxon_namespace, expected_taxon_labels):
                self.assertEqual(taxon.label, expected_label)

            # check general tree state
            self.assertEqual(num_trees, expected_num_trees)
            self.assertIs(is_rooted, expected_is_rooted)

            splits_ref = paupsplitsreference.get_splits_reference(
                    splits_filename=splits_filename,
                    key_column_index=0,
                    )
            self.assertEqual(len(splits_ref), len(bipartition_counts))
            self.assertEqual(len(splits_ref), len(bipartition_freqs))
            if is_rooted:
                splits_ref_bitmasks = set([splits_ref[x]["unnormalized_split_bitmask"] for x in splits_ref])
            else:
                splits_ref_bitmasks = set([splits_ref[x]["normalized_split_bitmask"] for x in splits_ref])
            counts_keys = set(bipartition_counts.keys())
            freqs_keys = set(bipartition_freqs.keys())
            self.assertEqual(len(counts_keys), len(splits_ref_bitmasks))
            self.assertEqual(counts_keys, splits_ref_bitmasks, "\n    {}\n\n    {}\n\n".format(sorted(counts_keys), sorted(splits_ref_bitmasks)))
            for split_str_rep in splits_ref:
                ref = splits_ref[split_str_rep]
                self.assertEqual(split_str_rep, ref["bipartition_string"])
                self.assertEqual(paup.PaupService.bipartition_groups_to_split_bitmask(split_str_rep, normalized=False),
                        ref["unnormalized_split_bitmask"])
                self.assertEqual(paup.PaupService.bipartition_groups_to_split_bitmask(split_str_rep, normalized=True),
                        ref["normalized_split_bitmask"])
                split_bitmask = paup.PaupService.bipartition_groups_to_split_bitmask(split_str_rep, normalized=not is_rooted)
                self.assertEqual(bipartition_counts[split_bitmask], ref["count"])
                # self.assertAlmostEqual(bipartition_freqs[split_bitmask], ref["frequency"])
                self.assertAlmostEqual(bipartition_freqs[split_bitmask], ref["frequency"], 2) # PAUP* 4.10b: no very precise

        def test_group1(self):
            cetacean_taxon_labels = [
                "Bos taurus",
                "Balaena mysticetus",
                "Balaenoptera physalus",
                "Cephalorhynchus eutropia",
                "Delphinapterus leucas",
                "Delphinus delphis",
                "Eschrichtius robustus",
                "Globicephala melas",
                "Inia geoffrensis",
                "Kogia breviceps",
                "Kogia simus",
                "Lagenorhynchus albirostris",
                "Lagenorhynchus obscurus",
                "Lissodelphis peronii",
                "Megaptera novaeangliae",
                "Mesoplodon europaeus",
                "Mesoplodon peruvianus",
                "Phocoena phocoena",
                "Phocoena spinipinnis",
                "Physeter catodon",
                "Tursiops truncatus",
                "Ziphius cavirostris",
            ]
            issue_mth_taxon_labels = ["T{:02d}".format(i) for i in range(1, 60)]
            sources = [
                    ("cetaceans.mb.no-clock.mcmc.trees"    , 251, False, False), # Trees explicitly unrooted
                    ("cetaceans.mb.no-clock.mcmc.weighted-01.trees" , 251, False , True), # Weighted
                    ("cetaceans.mb.no-clock.mcmc.weighted-02.trees" , 251, False , True), # Weighted
                    ("cetaceans.mb.no-clock.mcmc.weighted-03.trees" , 251, False , True), # Weighted
                    ("cetaceans.mb.strict-clock.mcmc.trees", 251, True , False), # Trees explicitly rooted
                    ("cetaceans.mb.strict-clock.mcmc.weighted-01.trees" , 251, True , True), # Weighted
                    ("cetaceans.mb.strict-clock.mcmc.weighted-02.trees" , 251, True , True), # Weighted
                    ("cetaceans.mb.strict-clock.mcmc.weighted-03.trees" , 251, True , True), # Weighted
                    ("cetaceans.raxml.bootstraps.trees"    , 250, True , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                    ("cetaceans.raxml.bootstraps.weighted-01.trees"    , 250, True , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                    ("cetaceans.raxml.bootstraps.weighted-02.trees"    , 250, True , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                    ("cetaceans.raxml.bootstraps.weighted-03.trees"    , 250, True , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
                    ("issue_mth_2009-02-03.rooted.nexus"   , 100, True , False), # 100 trees (frequency column not reported by PAUP)
                    ("issue_mth_2009-02-03.unrooted.nexus" , 100, False , False), # 100 trees (frequency column not reported by PAUP)
            ]
            splits_filename_template = "{stemname}.is-rooted-{is_rooted}.use-tree-weights-{use_weights}.burnin-{burnin}.splits.txt"
            for tree_filename, num_trees, treefile_is_rooted, treefile_is_weighted in sources:
                stemname = tree_filename
                if "cetacean" in tree_filename:
                    expected_taxon_labels = cetacean_taxon_labels
                    taxa_definition_filepath = pathmap.tree_source_path("cetaceans.taxa.nex")
                else:
                    expected_taxon_labels = issue_mth_taxon_labels
                    taxa_definition_filepath = pathmap.tree_source_path("issue_mth_2009-02-03.unrooted.nexus")
                for use_weights in (False, True, None):
                    for paup_read_as_rooted in (None, True, False):
                        for paup_burnin in (0, 150):
                            if tree_filename.startswith("issue_mth") and paup_burnin > 0:
                                continue
                            if paup_read_as_rooted is None:
                                expected_is_rooted = treefile_is_rooted
                            elif paup_read_as_rooted:
                                expected_is_rooted = True
                            else:
                                expected_is_rooted = False
                            splits_filename = splits_filename_template.format(
                                    stemname=stemname,
                                    is_rooted=paup_read_as_rooted,
                                    use_weights=use_weights,
                                    burnin=paup_burnin)
                            self.check_splits_counting(
                                    tree_filename=tree_filename,
                                    taxa_definition_filepath=taxa_definition_filepath,
                                    splits_filename=splits_filename,
                                    paup_as_rooted=paup_read_as_rooted,
                                    paup_use_tree_weights=use_weights,
                                    paup_burnin=paup_burnin,
                                    expected_taxon_labels=expected_taxon_labels,
                                    expected_is_rooted=expected_is_rooted,
                                    expected_num_trees=num_trees-paup_burnin)

    class PaupWrapperTaxaParse(ExtendedTestCase):

        def setUp(self):
            self.taxa_filepath = None
            self.expected_taxlabels = None

        def check_labels(self):
            p = paup.PaupRunner()
            p.stage_execute_file(self.taxa_filepath)
            p.stage_list_taxa()
            p.run()
            taxon_namespace = p.parse_taxon_namespace()
            self.assertEqual(len(taxon_namespace), len(self.expected_taxlabels))
            for i, t in enumerate(taxon_namespace):
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

#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Test the PAUP* wrapper
"""

import os
import sys
import csv
import collections
import unittest

from dendropy.test.support import pathmap
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.utility import container
from dendropy.utility import messaging
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
_LOG = messaging.get_logger(__name__)

from dendropy.calculate import treesplit
from dendropy.interop import paup

def get_splits_reference(
        splits_filename,
        splits_dir=None,
        key_column_index=0):
    # Key columns are:
    #     0   : PAUP* bipartition string representation '....**...' etc.
    #     1   : unnormalized split bitmask (for rooted trees)
    #     2   : normalized split bitmask (for unrooted trees)
    #     3   : (weighted) counts
    #     4   : (weighted) frequencies
    if splits_dir is not None:
        splits_filepath = os.path.join(splits_dir, splits_filename)
    else:
        splits_filepath = pathmap.splits_source_path(splits_filename)
    d = collections.OrderedDict()
    with open(splits_filepath, "r") as src:
        for row in src:
            content = row.split("#")[0]
            if not content:
                continue
            fields = content.split("\t")
            assert len(fields) == 5, "{}: {}".format(content, fields)
            key = fields[key_column_index]
            d[key] = {
                "bipartition_string": fields[0],
                "unnormalized_split_bitmask": fields[1],
                "normalized_split_bitmask": fields[2],
                "count": fields[3],
                "frequency": fields[4],
            }
    return d


if not paup.DENDROPY_PAUP_INTEROPERABILITY:
    _LOG.warn("PAUP interoperability not available: skipping PAUP tests")
else:

    class PaupWrapperRepToSplitMaskTest(unittest.TestCase):

        def setUp(self):
            self.ps = paup.PaupService()

        def testUnnormalized(self):
            for i in range(0xFF):
                s = treesplit.split_as_string(i, 8, ".", "*")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=False)
                self.assertEqual(r, i, "%s  =>  %s  =>  %s" \
                    % (treesplit.split_as_string(i, 8), s, treesplit.split_as_string(r, 8)))

        def testNormalized0(self):
            for i in range(0xFF):
                s = treesplit.split_as_string(i, 8, "*", ".")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=True)
                normalized = container.NormalizedBitmaskDict.normalize(i, 0xFF, 1)
                self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                    % (treesplit.split_as_string(i, 8), s, treesplit.split_as_string(normalized, 8)))

        def testNormalized1(self):
            for i in range(0xFF):
                s = treesplit.split_as_string(i, 8, ".", "*")[::-1]
                r = paup.PaupService.bipartition_groups_to_split_bitmask(s, normalized=True)
                normalized = container.NormalizedBitmaskDict.normalize(i, 0xFF, 1)
                self.assertEqual(r, normalized, "%s  =>  %s  =>  %s" \
                    % (treesplit.split_as_string(i, 8), s, treesplit.split_as_string(normalized, 8)))

    class PaupWrapperSplitsParse(ExtendedTestCase):

        def check_splits_counting(self,
                tree_filename,
                taxa_definition_filepath,
                splits_filename,
                paup_as_rooted,
                paup_ignore_tree_weights,
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
                    ignore_tree_weights=paup_ignore_tree_weights,
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

            splits_ref = get_splits_reference(
                    splits_filename=splits_filename,
                    key_column_index=0,
                    )
            self.assertEqual(len(splits_ref), len(bipartition_counts))
            self.assertEqual(len(splits_ref), len(bipartition_freqs))
            # for split_str_rep in self.expected_split_freqs:
            #     split_bitmask = paup.PaupService.bipartition_groups_to_split_bitmask(split_str_rep, normalized=not is_rooted)
            #     self.assertIn(split_bitmask, bipartition_counts)
            #     self.assertEqual(self.expected_split_freqs[split_str_rep], bipartition_counts[split_bitmask])
            # sd = paup.build_split_distribution(bipartition_counts,
            #                                    tree_count,
            #                                    taxon_namespace,
            #                                    is_rooted=is_rooted)
            # sf = sd.split_frequencies
            # for g in bipartition_counts:
            #     s = paup.paup_group_to_mask(g, normalized=not is_rooted)
            #     self.assertIn(s, sd.split_counts)
            #     self.assertEqual(sd.split_counts[s], bipartition_counts[g])
            #     self.assertEqual(sd.total_trees_counted, self.expected_num_trees)
            #     self.assertAlmostEqual(sf[s], float(bipartition_counts[g]) / self.expected_num_trees)


        def test_group1(self):
            expected_taxon_labels = [
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
            sources = [
                    ("cetaceans.mb.no-clock.mcmc.trees"    , 251, False, False), # Trees explicitly unrooted
                    ("cetaceans.mb.strict-clock.mcmc.trees", 251, True , False), # Trees explicitly rooted
                    ("cetaceans.raxml.bootstraps.trees"    , 250, True , False), # No tree rooting statement; PAUP defaults to rooted, DendroPy defaults to unrooted
            ]
            taxa_definition_filepath = pathmap.tree_source_path("cetaceans.taxa.nex")
            splits_filename_template = "{stemname}.trees.is-rooted-{is_rooted}.ignore-tree-weights-{ignore_weights}.burnin-{burnin}.splits.txt"
            for tree_filename, num_trees, treefile_is_rooted, treefile_is_weighted in sources:
                stemname = tree_filename.rsplit(".", 1)[0]
                for paup_read_as_rooted in (None, True, False):
                    for paup_burnin in (0, 150):
                        if paup_read_as_rooted is None:
                            expected_is_rooted = treefile_is_rooted
                        elif paup_read_as_rooted:
                            expected_is_rooted = True
                        else:
                            expected_is_rooted = False
                        splits_filename = splits_filename_template.format(
                                stemname=stemname,
                                is_rooted=paup_read_as_rooted,
                                ignore_weights=None,
                                burnin=paup_burnin)
                        self.check_splits_counting(
                                tree_filename=tree_filename,
                                taxa_definition_filepath=taxa_definition_filepath,
                                splits_filename=splits_filename,
                                paup_as_rooted=paup_read_as_rooted,
                                paup_ignore_tree_weights=None,
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

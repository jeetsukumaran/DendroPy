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
Tests of split distribution counting.
"""

import unittest
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import runlevel
from dendropy.test.support.extendedtest import ExtendedTestCase
from dendropy.utility import messaging
from dendropy.interop import paup
from dendropy import dataio
from dendropy import treesum
from dendropy import treesplit
import dendropy

_LOG = messaging.get_logger(__name__)

if not paup.DENDROPY_PAUP_INTEROPERABILITY:
    _LOG.warn("PAUP interoperability not available: skipping split counting tests")
else:

    class SplitCountTest(ExtendedTestCase):

        def setUp(self):
            self.test_cases = [('pythonidae.reference-trees.nexus', 'pythonidae.reference-trees.nexus')]
            if runlevel.is_test_enabled(runlevel.SLOW, _LOG, self.__class__.__name__):
                self.test_cases.extend([
                    ('feb032009.trees.nexus', 'feb032009.trees.nexus'),
                    ('maj-rule-bug1.trees.nexus', 'maj-rule-bug1.trees.nexus'),
                    ('maj-rule-bug2.trees.nexus', 'maj-rule-bug2.trees.nexus'),
                ])

        def countSplits(self, tc, is_rooted):
            _LOG.info(tc[0] + "; " + tc[1])
            tree_filepaths = [pathmap.tree_source_path(tc[0])]
            taxa_filepath = pathmap.tree_source_path(tc[1])
            paup_sd = paup.get_split_distribution(
                    tree_filepaths,
                    taxa_filepath,
                    is_rooted=is_rooted,
                    burnin=0)
            taxon_set = paup_sd.taxon_set
            dp_sd = treesplit.SplitDistribution(taxon_set=taxon_set)
            dp_sd.ignore_edge_lengths = True
            dp_sd.ignore_node_ages = True
            dp_sd.is_rooted = is_rooted

            _LOG.debug("Taxon set: %s" % [t.label for t in taxon_set])
            taxa_mask = taxon_set.all_taxa_bitmask()
            taxon_set.lock()
            for tree_filepath in tree_filepaths:
                for tree in dataio.tree_source_iter(
                        stream=open(tree_filepath, "rU"),
                        schema='nexus',
                        taxon_set=taxon_set,
                        as_rooted=is_rooted):
                    self.assertIs(tree.taxon_set, dp_sd.taxon_set)
                    self.assertIs(tree.taxon_set, taxon_set)
                    treesplit.encode_splits(tree)
                    dp_sd.count_splits_on_tree(tree)

            self.assertEqual(dp_sd.total_trees_counted, paup_sd.total_trees_counted)

            # SplitsDistribution counts trivial splits, whereas PAUP*
            # contree does not, so the following will not work
    #            assert len(dp_sd.splits) == len(paup_sd.splits),\
    #                 "dp = %d, sd = %d" % (len(dp_sd.splits), len(paup_sd.splits))

            taxa_mask = taxon_set.all_taxa_bitmask()
            for split in dp_sd.splits:
                if not treesplit.is_trivial_split(split, taxa_mask):
                    self.assertIn(split, paup_sd.splits)
                    self.assertEqual(dp_sd.split_counts[split], paup_sd.split_counts[split])
                    paup_sd.splits.remove(split)

            # if any splits remain, they were not
            # in dp_sd or were trivial
            remaining_splits = list(paup_sd.splits)
            for split in remaining_splits:
                if treesplit.is_trivial_split(split, taxa_mask):
                    paup_sd.splits.remove(split)
            self.assertEqual(len(paup_sd.splits), 0)

        def testUnrootedSplitCounts(self):
            for tc in self.test_cases:
                self.countSplits(tc, is_rooted=False)

        def testRootedSplitCounts(self):
            for tc in self.test_cases:
                self.countSplits(tc, is_rooted=True)

    if __name__ == "__main__":
        unittest.main()

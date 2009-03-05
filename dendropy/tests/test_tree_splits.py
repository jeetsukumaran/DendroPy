#! /usr/bin/env python

############################################################################
##  test_tree_splits.py
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
Test split summarization, consensus tree building etc.
"""

import unittest
import subprocess
import re
import os
import sys

from dendropy import nexus
from dendropy import get_logger
import dendropy.tests
from dendropy.tests import paup, is_test_enabled, TestLevel

_LOG = get_logger("Splits")

from dendropy import nexus
from dendropy import treesum
from dendropy import taxa
from dendropy import splits

large_cases = [ #('7180.tre', '7180.tre'),
                     #('terrarana.random.unrooted.100.tre', 'terrarana.random.unrooted.100.tre'),
                     ('terrarana.random.unrooted.30.tre', 'terrarana.random.rooted.30.tre'),
                     ('anolis.mcmct.trees.nexus', 'anolis.chars.nexus'),
]
small_cases = [ ('feb032009.tre', 'feb032009.tre'),
                     ('maj-rule-bug1.tre', 'maj-rule-bug1.tre'),
                     ('maj-rule-bug2.tre', 'maj-rule-bug2.tre'),
]
med_cases = [ ('primates.mcmct.trees.nexus', 'primates.chars.nexus'),
]
test_cases = small_cases
if is_test_enabled(TestLevel.NORMAL, _LOG, module_name=__name__, message="skipping medium tree files"):
    test_cases.extend(med_cases)
if is_test_enabled(TestLevel.SLOW, _LOG, module_name=__name__, message="skipping large tree files"):
    test_cases.extend(large_cases)

class SplitFreqsTest(unittest.TestCase):


    def testSplits(self):
        unrooted = True
        for tc in test_cases:
            tree_filepaths = [dendropy.tests.data_source_path(tc[0])]
            taxa_filepath = dendropy.tests.data_source_path(tc[1])
            paup_sd = paup.get_split_distribution(tree_filepaths, taxa_filepath,
                        unrooted=unrooted, burnin=0)
            taxa_block = paup_sd.taxa_block
            dp_sd = splits.SplitDistribution(taxa_block=taxa_block)
            dp_sd.ignore_edge_lengths = True
            dp_sd.ignore_node_ages = True
            dp_sd.unrooted = unrooted

            taxa_mask = taxa_block.all_taxa_bitmask()
            taxa_block.lock()
            for tree_filepath in tree_filepaths:
                for tree in nexus.iterate_over_trees(open(tree_filepath, "rU"), taxa_block=taxa_block):
                    #_LOG.debug("tree = %s" % str(tree))
                    splits.encode_splits(tree)
                    dp_sd.count_splits_on_tree(tree)

            self.assertEqual(dp_sd.total_trees_counted, paup_sd.total_trees_counted)

            # SplitsDistribution counts trivial splits, whereas PAUP*
            # contree does not, so the following will not work
#             assert len(dp_sd.splits) == len(paup_sd.splits),\
#                 "dp = %d, sd = %d" % (len(dp_sd.splits), len(paup_sd.splits))

            taxa_mask = taxa_block.all_taxa_bitmask()
            for split in dp_sd.splits:
                if not splits.is_trivial_split(split, taxa_mask):
                    self.assertTrue(split in paup_sd.splits)
                    self.assertEqual(dp_sd.split_counts[split], paup_sd.split_counts[split])
                    paup_sd.splits.remove(split)

            # if any splits remain here, they were not
            # in dp_sd
            assert len(paup_sd.splits) == 0


if __name__ == "__main__":
    unittest.main()

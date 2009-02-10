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
from dendropy.tests import paup

_LOG = get_logger("Splits")


from dendropy import treesum
from dendropy import taxa   
from dendropy import splits

class SplitsOnTreesTest(unittest.TestCase):
    
    def setUp(self):
        self.data_files = [   
#             ["feb032009.tre", "feb032009.tre"],        
            ["anolis.chars.nexus", "anolis.mcmct.trees.nexus"],
            ["primates.chars.nexus", "primates.mcmct.trees.newick"],
        ]
                
    def check_split_summarization(self, data_file, tree_file, min_clade_freq=0.5, burnin=0):
    
        # get PAUP's version of the splits ...
        _LOG.info("PAUP is counting ...")
        nexus_tree_file = tree_file.replace("newick", "nexus")
        tax_labels, paup_biparts, paup_biparts_count, paup_biparts_freqs = paup.bipartitions(data_file, nexus_tree_file, min_clade_freq, burnin)
        paup_biparts_complemented = []
        for bipart in paup_biparts:
            paup_biparts_complemented.append(bipart.replace('.','X').replace('*','.').replace('X','*'))
        
        # get our version ...
        _LOG.info("DendroPy is counting ...")
        taxa_block = taxa.TaxaBlock() 
        for tax_label in tax_labels:
            taxa_block.add_taxon(label=tax_label.replace(' ', '_'))                  
        sd = splits.SplitDistribution(taxa_block=taxa_block)        
        tsum = treesum.TreeSummarizer()
        tsum.verbose = True
        tsum.write_message = _LOG.debug
        tsum.progress_message_suffix = " (%s)" % os.path.basename(tree_file)
        sd = tsum.count_splits_on_trees([tree_file], 
                               nexus.iterate_over_trees, 
                               split_distribution=sd)  
        contree = tsum.tree_from_splits(sd)
        
        
        _LOG.info("Collating information ...")
        dendropy_split_strings = []
        dendropy_string_split_map = {}
        dendropy_complemented_split_strings = []
        dendropy_string_complemented_split_map = {}
        taxa_mask = taxa_block.all_taxa_bitmask()
        for split in sd.splits:
            if splits.is_non_singleton_split(split, taxa_mask):
                rsplit = taxa_block.complement_split_bitmask(split) 
                if rsplit:
                    ss = splits.split_as_string_rev(split, len(sd.taxa_block), '.', '*')
                    dendropy_split_strings.append(ss)
                    dendropy_string_split_map[ss] = split
                    ss = splits.split_as_string_rev(rsplit, len(sd.taxa_block), '.', '*')
                    dendropy_complemented_split_strings.append(ss)
                    dendropy_string_complemented_split_map[ss] = rsplit
                
        # make sure the distinct splits are the same across both versions
        _LOG.info("Checking for correspondence in split identity ...")
        for s in dendropy_split_strings:
            assert (s in paup_biparts) or (s in paup_biparts_complemented), \
                            "PAUP did not find: %s" % s
        for s in paup_biparts:
            assert (s in dendropy_split_strings) \
                    or (s in dendropy_complemented_split_strings), \
                            "DendroPy did not find: %s" % s
                            
        # make sure the counts/freqs are the same
        _LOG.info("Checking for correspondence in split counts/frequences ...")        
        split_freqs = sd.calc_freqs()
        for idx, s in enumerate(paup_biparts):
            if s in dendropy_split_strings:
                split = dendropy_string_split_map[s]
                count = sd.split_counts[split]
                freq = split_freqs[split]
            else:
                s2 = paup_biparts_complemented[idx]
                split = dendropy_string_split_map[s2]
                count = sd.split_counts[split]
                freq = split_freqs[split]
                
            _LOG.debug("PAUP: %d (%f), DendroPy: %d (%f)" % (paup_biparts_count[s],
                                                             paup_biparts_freqs[s],
                                                             count,
                                                             freq*100))
            
            assert paup_biparts_count[s] == count, \
                   "Counts of split '%s': Expecting %d, but found %d" \
                    % (s, paup_biparts_count[s], count)
            assert dendropy.tests.is_almost_equal(paup_biparts_freqs[s], freq), \
                   "Frequency of split '%s': Expecting %f, but found %f" \
                    % (s, paup_biparts_freqs[s], perc)
                    
        _LOG.info("\n--SUCCESS--\n")                             
    def testSplitsSummary(self):
        for df in self.data_files:
            char_file = dendropy.tests.data_source_path(df[0])
            tree_file = dendropy.tests.data_source_path(df[1])
            if os.path.exists(char_file) and os.path.exists(tree_file):
                _LOG.info("Checking splits in: %s" % os.path.basename(tree_file))
                self.check_split_summarization(data_file=char_file, tree_file=tree_file)
                
if __name__ == "__main__":
    unittest.main()

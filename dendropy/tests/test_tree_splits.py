#! /usr/bin/env python

############################################################################
##  test_tree_splits.py
##
##  Part of the DendroPy phylogenetic computation library.
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
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
            ["anolis.chars.nexus", "anolis.mbcon.trees.nexus"],
            ["anolis.chars.nexus", "anolis.mcmct.trees.nexus"],
            ["primates.chars.nexus", "primates.mcmct.trees.nexus"],
        ]
                
    def check_split_summarization(self, data_file, tree_file, min_clade_freq=0.5, burnin=0):
        tax_labels, biparts, bpc, bpf = paup.bipartitions(data_file, tree_file, min_clade_freq, burnin)
        biparts = bpc.keys()
        biparts_c = []
        for bipart in biparts:
            biparts_c.append(bipart.replace('.','X').replace('*','.').replace('X','*'))
        
        # create taxa
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
        dendropy_split_strings = []
        dendropy_split_strings_c = []
        for split in sd.splits:
            if splits.is_non_singleton_split(split) and (split ^ taxa_block.all_taxa_bitmask()):
                dendropy_split_strings.append(splits.split_as_string_rev(split, sd.taxa_block, '.', '*'))
                dendropy_split_strings_c.append(splits.split_as_string_rev(split, sd.taxa_block, '*', '.'))
      
        for s in dendropy_split_strings:
            self.failUnless(s in biparts or s in biparts_c,
                            "PAUP did not find: %s" % s)       
        for s in biparts:
            self.failUnless(s in dendropy_split_strings or s in dendropy_split_strings_c,
                            "DendroPy did not find: %s" % s)
                            
        _LOG.info("\n--SUCCESS--\n")                            
                            
    def test_splits_summary(self):
        for df in self.data_files:
            char_file = dendropy.tests.data_source_path(df[0])
            tree_file = dendropy.tests.data_source_path(df[1])
            if os.path.exists(char_file) and os.path.exists(tree_file):
                _LOG.info("Checking splits in: %s" % os.path.basename(tree_file))
                self.check_split_summarization(data_file=char_file, tree_file=tree_file)
                
if __name__ == "__main__":
    unittest.main()

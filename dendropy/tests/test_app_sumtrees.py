#! /usr/bin/env python

############################################################################
##  test_app_sumtrees.py
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
Put SumTrees through its paces.
"""

import unittest
import subprocess
import tempfile
import re
import os
import sys

from dendropy import nexus
from dendropy import get_logger
import dendropy.tests
from dendropy.tests import is_test_enabled, TestLevel
import itertools
from dendropy.treedists import symmetric_difference
from dendropy.splits import SplitDistribution
from dendropy.treesum import TreeSummarizer
from dendropy import dataio
from dendropy.nexus import RootingInterpretation
from dendropy.splits import encode_splits
from dendropy.datasets import Dataset
_LOG = get_logger("SumTreesTesting")


class SumTreesTest(unittest.TestCase):
    
    def setUp(self):                                   
        self.sumtrees_path = dendropy.tests.scripts_source_path("sumtrees.py")
                
    def compose_sumtrees_command(self, args):
        return "%s %s" % (self.sumtrees_path, " ".join(args))        
        
    def runSumTrees(self, args):
        command = self.compose_sumtrees_command(args)
        _LOG.debug("Invocation:\n"+command)
        run = subprocess.Popen(command, 
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        stdout, stderr = run.communicate()
        _LOG.debug(stderr)
        _LOG.debug(stdout)
        assert run.returncode == 0, "\nSumTrees exited with error: %s" % stderr

    def test3Feb2009MajRuleBug(self):
        if not is_test_enabled(TestLevel.NORMAL, _LOG, module_name=__name__, message="skipping sumtree argument tests"):
            return
        fn1 = dendropy.tests.data_source_path("maj-rule-bug1.tre")
        fn2 = dendropy.tests.data_source_path("maj-rule-bug2.tre")
        d = Dataset()
        
        tb1 = d.read_trees(open(fn1,"rU"), format="NEXUS", encode_splits=True, rooted=RootingInterpretation.UNROOTED)
        tb2 = d.read_trees(open(fn2,"rU"), format="NEXUS", encode_splits=True, rooted=RootingInterpretation.UNROOTED)
        taxa1 = d.taxa_blocks[0]
        self.assertEqual(taxa1, tb2[0].taxa_block)
        
        firstSD = SplitDistribution(taxa_block=taxa1)
        secondSD = SplitDistribution(taxa_block=taxa1)
        
        for o,t in itertools.izip(tb1, tb2):
            #encode_splits(o)
            #encode_splits(t)
            firstSD.count_splits_on_tree(o)
            secondSD.count_splits_on_tree(t)

        ts = TreeSummarizer()
        n_times = 1 # keep set to 1 except for benchmarking tree_from_splits
        for i in xrange(n_times):
            firstMR = ts.tree_from_splits(firstSD, min_freq=0.5)
            secondMR = ts.tree_from_splits(secondSD, min_freq=0.5)
        self.assertEqual(0, symmetric_difference(firstMR, secondMR))

    def testSumTreeOptions(self):
        support_file = dendropy.tests.data_source_path("anolis.mbcon.trees.nexus")
        target_file = dendropy.tests.data_source_path("anolis.mbcon.trees.nexus")
        outfile = tempfile.mktemp()
        
        # default options,
        self.runSumTrees([support_file])
        if not is_test_enabled(TestLevel.SLOW, _LOG, module_name=__name__, message="skipping sumtree argument tests"):
            return
        # default options, multiple files
        self.runSumTrees([support_file, support_file, support_file])
                       
        # burnin 
        self.runSumTrees(["--burnin=100", support_file])
        
        # target tree
        self.runSumTrees(["--target=%s" % target_file, support_file])
        
        # 95% consensus tree
        self.runSumTrees(["--min-clade-freq=0.95", support_file])
        
        # no branch lengths
        self.runSumTrees(["--no-branch-lengths", support_file])
        
        # support as labels
        self.runSumTrees(["--support-as-labels", support_file])        
        
        # support as branch lengths
        self.runSumTrees(["--support-as-lengths", support_file])
        
        # support as percentages
        self.runSumTrees(["--percentages", support_file])
        
        # support decimals
        self.runSumTrees(["--decimals=0", support_file])        
        
        # output to tmp file
        if os.path.exists(outfile):
            os.remove(outfile)
        self.runSumTrees(["--output=%s" % outfile, support_file])
        self.runSumTrees(["--replace --output=%s" % outfile, support_file])
        if os.path.exists(outfile):
            os.remove(outfile)        
        
        # no taxa block
        self.runSumTrees(["--no-taxa-block", support_file])
        
        # no taxa block
        self.runSumTrees(["--no-meta-comments", support_file])
        
        # additional comments
        self.runSumTrees(["-m 'Test run of SumTrees'", support_file])
        
        # newick format
        self.runSumTrees(["--to-newick", support_file])
        
        # ignore missing support files
        self.runSumTrees(["--ignore-missing-support", support_file, "dummy"])
        
        # ignore missing target
        self.runSumTrees(["--ignore-missing-target", "--target=dummy", support_file])
        
        # quiet
        self.runSumTrees(["--quiet", support_file])        
        
        
if __name__ == "__main__":
    unittest.main()

#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Tests SumTrees.
"""

import os
import unittest
import subprocess
from dendropy.test.support import pathmap
from dendropy.test.support import runlevel
from dendropy.utility.messaging import get_logger
from dendropy.test.support import extendedtest

import dendropy

_LOG = get_logger(__name__)

class SumTreesTester(extendedtest.ExtendedTestCase):

    def setUp(self, sumtrees_path=None):
        if sumtrees_path is None:
            self.sumtrees_path = pathmap.script_source_path(os.path.join('sumtrees', 'sumtrees.py'))
        else:
            self.sumtrees_path = sumtrees_path

    def execute_sumtrees(self, args):
        if isinstance(args, str):
            args = [args]
        cmd = [self.sumtrees_path] + args
        _LOG.debug(" ".join(cmd))
        p = subprocess.Popen(cmd,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        retcode = p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()
        return retcode, stdout, stderr

    def assertExecuteWithoutError(self, cmd, clean_stderr=True):
        retcode, stdout, stderr = self.execute_sumtrees(cmd)
        self.assertEqual(retcode, 0)
        if clean_stderr:
            self.assertEqual(stderr, "")

class SumTreesNodeAgesTester(SumTreesTester):


    def testSummarizeNodeAgesOnMCCT(self):
        """
        SumTrees: summarizing node ages on MCCT topology.
        """
        if runlevel.is_test_enabled(runlevel.EXHAUSTIVE, _LOG, self.__class__.__name__):
            path_to_src = pathmap.tree_source_path("primates.beast.mcmc.trees")
            path_to_target = pathmap.tree_source_path("primates.beast.mcct.noedgelens.tree")
            args = ["-b",
                    "40",
                    "-e",
                    "mean-age",
                    "-t",
                    path_to_target,
                    path_to_src]
            retcode, stdout, stderr = self.execute_sumtrees(args)
            self.assertEqual(retcode, 0)

            taxa = dendropy.TaxonSet()
            exp_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path("primates.beast.mcct.meanh.tre"), "nexus", taxon_set=taxa)
            obs_tree = dendropy.Tree.get_from_string(stdout, "nexus", taxon_set=taxa)
            exp_tree.update_splits()
            exp_tree.calc_node_ages()
            obs_tree.update_splits()
            obs_tree.calc_node_ages()
            self.assertEqual(exp_tree.split_edges.keys(), obs_tree.split_edges.keys())
            splits = exp_tree.split_edges.keys()
            for split in splits:
                exp_edge = exp_tree.split_edges[split]
                obs_edge = obs_tree.split_edges[split]
                self.assertAlmostEqual(obs_edge.head_node.age, exp_edge.head_node.age)
        else:
            _LOG.info("Skipping test (set 'DENDROPY_TESTING_LEVEL=EXHAUSTIVE' to run)")


class SumTreeHelpAndVersionTests(SumTreesTester):

    def testVersion(self):
        self.assertExecuteWithoutError("--version")

    def testHelp1(self):
        self.assertExecuteWithoutError("-h")

    def testHelp2(self):
        self.assertExecuteWithoutError("--help")

if __name__ == "__main__":
    unittest.main()

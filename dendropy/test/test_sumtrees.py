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
from dendropy.utility.messaging import get_logger

_LOG = get_logger(__name__)

class SumTreesTester(unittest.TestCase):

    def setUp(self, sumtrees_path=None):
        if sumtrees_path is None:
            self.sumtrees_path = pathmap.script_source_path(os.path.join('sumtrees', 'sumtrees.py'))
        else:
            self.sumtrees_path = sumtrees_path

    def execute_sumtrees(self, args):
        if isinstance(args, str):
            args = [args]
        cmd = [self.sumtrees_path] + args
        _LOG.info(cmd)
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

class SumTreeHelpAndVersionTests(SumTreesTester):

    def testVersion(self):
        self.assertExecuteWithoutError("--version")

    def testHelp1(self):
        self.assertExecuteWithoutError("-h")

    def testHelp2(self):
        self.assertExecuteWithoutError("--help")

if __name__ == "__main__":
    unittest.main()

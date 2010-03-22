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

class SumTreesOptionsTests(SumTreesTester):

    def testHelp(self):
        self.assertExecuteWithoutError("-h")


if __name__ == "__main__":
    unittest.main()

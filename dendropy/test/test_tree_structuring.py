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
Tests native tree structuring routines.
"""

import unittest
from dendropy.test.support import pathmap
from dendropy.utility import messaging
import dendropy

_LOG = messaging.get_logger(__name__)

class TestTreeRestructures(unittest.TestCase):

    def testLadderizeLeft(self):
        """
        Dummy test!
        """
        tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        tree.ladderize()

    def testLadderizeRight(self):
        """
        Dummy test!
        """
        tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        tree.ladderize(right=True)

if __name__ == "__main__":
    unittest.main()


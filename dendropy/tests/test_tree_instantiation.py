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
Tests composition and traversal of trees.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import messaging
import dendropy.tests
from dendropy import Tree, TaxonSet

_LOG = messaging.get_logger(__name__)

class TreeInstantiationTest(unittest.TestCase):

    def test_tree_build(self):
        # from file, w/keywords
        t1 = Tree(istream=StringIO("((A,B),(C,D));"), format="newick")


if __name__ == "__main__":
    unittest.main()

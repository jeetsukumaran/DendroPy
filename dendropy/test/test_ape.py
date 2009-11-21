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
Tests the interaction with APE.
"""

import unittest
from dendropy.utility import error
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.interop import ape
import dendropy

class TreesRoundTrip(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.trees = datagen.reference_tree_list()

    def testApeTreeRoundTrip(self):
        for i, t1 in enumerate(self.trees):
            ape_t = ape.as_ape_tree(t1)
            rt = ape.as_dendropy_tree(ape_t, taxon_set=t1.taxon_set)
            self.assertTrue(isinstance(rt, dendropy.Tree), type(rt))
            self.assertDistinctButEqual(
                    t1, rt, distinct_taxa=False)

if __name__ == "__main__":
    unittest.main()


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
from dendropy.test.support import pathmap
from dendropy.interop import ete
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

import dendropy

if not ete.DENDROPY_ETE_INTEROPERABILITY:
    _LOG.warn("ETE interoperability not available: skipping ETE tests")
else:
    class DataRoundTrip(datatest.DataObjectVerificationTestCase):

        def setUp(self):
            self.trees = datagen.reference_tree_list()

        def testTreeRoundTrip(self):
            t1 = self.trees[0]
            ete_t = ete.as_ete_object(t1)
            rt = ete.as_dendropy_object(ete_t, taxon_set=t1.taxon_set)
            self.assertTrue(isinstance(rt, dendropy.Tree), type(rt))
            self.assertDistinctButEqual(t1, rt, distinct_taxa=False)

        def testTreeListRoundTrip(self):
            ete_t = ete.as_ete_object(self.trees)
            rt = ete.as_dendropy_object(ete_t, taxon_set=self.trees.taxon_set)
            self.assertTrue(isinstance(rt, dendropy.TreeList), type(rt))
            self.assertDistinctButEqual(self.trees, rt, distinct_taxa=False)

    if __name__ == "__main__":
        unittest.main()


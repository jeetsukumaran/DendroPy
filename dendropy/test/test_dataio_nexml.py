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
NEXUS data read/write parse/format tests.
"""

import sys
import os
import unittest
import tempfile

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
import dendropy
from dendropy.dataio import nexml
from dendropy import splitmask
from dendropy import treecalc


class NexmlRoundTripTest(datatest.DataObjectVerificationTestCase):

    ## Trees need special attention
    def assertDistinctButEqualTree(self, tree1, tree2, **kwargs):
        otaxa = tree1.taxon_set
        ts = dendropy.TaxonSet()
        tree1.reindex_taxa(ts, clear=True)
        tree2.reindex_taxa(ts)
        self.assertSame(tree1.taxon_set, tree2.taxon_set)
        self.assertNotSame(tree1.taxon_set, otaxa)
        self.assertDistinctButEqual(tree1.taxon_set, otaxa, **kwargs)
        splitmask.encode_splits(tree1)
        splitmask.encode_splits(tree2)
        rfdist = treecalc.robinson_foulds_distance(tree1, tree2)
        self.assertAlmostEqual(rfdist, 0)

    def testRoundTreeJustTrees(self):
        ds = dendropy.DataSet(datagen.reference_tree_list())
        self.roundTripDataSetTest(ds, "nexml", ignore_taxon_order=True)

    def testRoundTripReference(self):
        reference_dataset = datagen.reference_single_taxonset_dataset()
        self.roundTripDataSetTest(reference_dataset, "nexml", ignore_taxon_order=True)

    def testRoundTripProtein(self):
        s = pathmap.char_source_stream("caenophidia_mos.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTripStandard1(self):
        s = pathmap.char_source_stream("angiosperms.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTripStandard2(self):
        s = pathmap.char_source_stream("apternodus.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="nexus")
        for ca in d1.char_arrays:
            ca.markup_as_sequences = False
        self.roundTripDataSetTest(d1, "nexml")

if __name__ == "__main__":
    unittest.main()

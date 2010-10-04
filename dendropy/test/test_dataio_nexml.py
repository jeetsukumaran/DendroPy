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
from dendropy import treesplit
from dendropy import treecalc


class NexmlRoundTripTest(datatest.DataObjectVerificationTestCase):

    ## Trees need special attention
    def assertDistinctButEqualTree(self, tree1, tree2, **kwargs):
        otaxa = tree1.taxon_set
        ts = dendropy.TaxonSet()
        tree1.reindex_taxa(ts, clear=True)
        tree2.reindex_taxa(ts)
        self.assertIs(tree1.taxon_set, tree2.taxon_set)
        self.assertIsNot(tree1.taxon_set, otaxa)
        self.assertDistinctButEqual(tree1.taxon_set, otaxa, **kwargs)
        treesplit.encode_splits(tree1)
        treesplit.encode_splits(tree2)
        rfdist = treecalc.robinson_foulds_distance(tree1, tree2)
        self.assertAlmostEqual(rfdist, 0)

    def testRoundTripProtein(self):
        s = pathmap.char_source_stream("caenophidia_mos.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTreeJustTrees(self):
        ds = dendropy.DataSet(datagen.reference_tree_list())
        self.roundTripDataSetTest(ds, "nexml", ignore_taxon_order=True)

    def testRoundTripReference(self):
        reference_dataset = datagen.reference_single_taxonset_dataset()
        self.roundTripDataSetTest(reference_dataset, "nexml", ignore_taxon_order=True)


    def testRoundTripStandard1(self):
        s = pathmap.char_source_stream("angiosperms.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTripStandard2(self):
        s = pathmap.char_source_stream("apternodus.chars.nexus")
        d1 = dendropy.DataSet(stream=s, schema="nexus")
        for ca in d1.char_matrices:
            ca.markup_as_sequences = False
        self.roundTripDataSetTest(d1, "nexml")

if __name__ == "__main__":
    unittest.main()

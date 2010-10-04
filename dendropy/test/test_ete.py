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


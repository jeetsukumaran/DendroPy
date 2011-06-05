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
Tests the interaction with Biopython.
"""

import unittest
from dendropy.utility import error
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.test.support import pathmap
from dendropy.interop import biopython
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

import dendropy

if not biopython.DENDROPY_BIOPYTHON_INTEROPERABILITY:
    _LOG.warn("Biopython interoperability not available: skipping Biopython tests")
else:
    class DataRoundTrip(datatest.DataObjectVerificationTestCase):

        def setUp(self):
            self.taxon_set = datagen.reference_taxon_set()
            self.dna_matrix = datagen.reference_dna_matrix(self.taxon_set)
            self.std_matrix = datagen.reference_standard_matrix(self.taxon_set)
            self.cont_matrix = datagen.reference_continuous_matrix(self.taxon_set)

        def testUnsupportedAlphabet(self):
            self.assertRaises(ValueError, biopython.as_biopython_object, self.cont_matrix)
            self.assertRaises(ValueError, biopython.as_biopython_object, self.std_matrix)


    if __name__ == "__main__":
        unittest.main()


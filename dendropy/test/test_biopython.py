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


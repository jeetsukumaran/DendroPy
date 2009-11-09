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
Tests creation, reading, update, deletion of CharArray objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import datatest
from dendropy.test.support import datagen
import dendropy

class DnaArrayTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.char_array1 = datagen.reference_dna_array()

    def testFromDnaCharArray(self):
        ca2 = dendropy.DnaCharacterArray(self.char_array1)
        self.assertDistinctButEqual(
            self.char_array1,
            ca2,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=False,
            equal_oids=False)

class StandardArrayTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.char_array1 = datagen.reference_standard_array()

    def testFromStandardCharArray(self):
        ca2 = dendropy.StandardCharacterArray(self.char_array1)
        self.assertDistinctButEqual(
            self.char_array1,
            ca2,
            char_type=dendropy.StandardCharacterArray,
            distinct_state_alphabets=True,
            distinct_taxa=False,
            equal_oids=False)

if __name__ == "__main__":
    unittest.main()

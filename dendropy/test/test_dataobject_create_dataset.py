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
Tests creation, reading, update, deletion of DataSet objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import datatest
from dendropy.test.support import datagen
import dendropy

class DataSetCreateTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.dataset = datagen.reference_single_taxonset_dataset()

    def testFromReadNexus(self):
        ds_str = self.dataset.as_string(format="nexus")
        ds2 = dendropy.DataSet(stream=StringIO(ds_str), format="nexus")
        self.assertDistinctButEqual(self.dataset, ds2, ignore_columns=True)

    def testFromCopy(self):
        ds2 = dendropy.DataSet(self.dataset)
        self.assertDistinctButEqual(self.dataset, ds2)

    def testSimpleCopyDnaArray(self):
        char_array = datagen.reference_dna_array()
        ds1 = dendropy.DataSet(char_array)
        self.assertEqual(len(ds1.char_arrays), 1)
        self.assertSame(ds1.char_arrays[0], char_array)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

    def testSimpleCopyStandardArray(self):
        char_array = datagen.reference_standard_array()
        ds1 = dendropy.DataSet(char_array)
        self.assertEqual(len(ds1.char_arrays), 1)
        self.assertSame(ds1.char_arrays[0], char_array)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

if __name__ == "__main__":
    unittest.main()

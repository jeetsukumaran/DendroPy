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
Tests FASTA I/O.
"""

import sys
import tempfile
import unittest
from cStringIO import StringIO
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
import dendropy

class TestFasta(datatest.DataObjectVerificationTestCase):

    def testAsStrReading(self):
        dataset = dendropy.DataSet(
                stream=open(pathmap.char_source_path("bad_names.fasta"), "rU"),
                schema='fasta',
                data_type='dna',
                row_type='str'
        )
        taxon_set = dataset.taxon_sets[0]
        label = [i.label for i in taxon_set]
        expected = ['a Bad name', 'another', 'a Badn,ame', 'a  nothe++-_=+r', 'an!@#$o^&*()}{_ther']
        self.assertEquals(label, expected)

    def testReadingAndWritingDataSet(self):
        ds1 = dendropy.DataSet(datagen.reference_dna_matrix())
        dataset = self.roundTripDataSetTest(ds1, "fasta", reader_kwargs={'data_type': 'dna'})

    def testReadingAndWritingCharMatrix(self):
        dna1 = datagen.reference_dna_matrix()
        output_path = pathmap.named_output_path(filename="roundtrip_test.fasta", suffix_timestamp=True)
        dna1.write_to_path(output_path, 'fasta')
        dna2 = dendropy.DnaCharacterMatrix.get_from_path(output_path, 'fasta')
        self.assertDistinctButEqual(dna1, dna2)

if __name__ == "__main__":
    unittest.main()


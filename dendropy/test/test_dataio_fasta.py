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
                format='DNAFasta',
                row_type='str'
        )
        taxon_set = dataset.taxon_sets[0]
        label = [i.label for i in taxon_set]
        expected = ['a Bad name', 'another', 'a Badn,ame', 'a  nothe++-_=+r', 'an!@#$o^&*()}{_ther']
        self.assertEquals(label, expected)

    def testReadingAndWriting(self):
        ds1 = dendropy.DataSet(datagen.reference_dna_array())
        dataset = self.roundTripDataSetTest(ds1, "DNAFASTA")

    def testRoundTripProtein(self):
        s = pathmap.char_source_stream("caenophidia_mos.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="fasta")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTripStandard1(self):
        s = pathmap.char_source_stream("angiosperms.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="fasta")
        self.roundTripDataSetTest(d1, "nexml")

    def testRoundTripStandard2(self):
        s = pathmap.char_source_stream("apternodus.chars.nexus")
        d1 = dendropy.DataSet(stream=s, format="fasta")
        for ca in d1.char_arrays:
            ca.markup_as_sequences = False
        self.roundTripDataSetTest(d1, "nexml")

if __name__ == "__main__":
    unittest.main()


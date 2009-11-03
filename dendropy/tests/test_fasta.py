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
from dendropy.tests import data_source_path
import dendropy

class TestFasta(unittest.TestCase):

    def testAsStrReading(self):
        dataset = dendropy.Dataset(
                istream=open(data_source_path("bad_names.fasta"), "rU"),
                format='DNAFasta',
                row_type='str'
        )
        taxon_set = dataset.taxon_sets[0]
        label = [i.label for i in taxon_set]
        expected = ['a Bad name', 'another', 'a Badn,ame', 'a  nothe++-_=+r', 'an!@#$o^&*()}{_ther']
        self.assertEquals(label, expected)

    def testAsStrReadingAndWriting(self):
        dataset = dendropy.Dataset(
                istream=open(data_source_path("bad_names.fasta"), "rU"),
                format="DNAFasta",
                row_type='str'
        )
        op = tempfile.TemporaryFile()
        dataset.write(ostream=op, format="FASTA")

if __name__ == "__main__":
    unittest.main()


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


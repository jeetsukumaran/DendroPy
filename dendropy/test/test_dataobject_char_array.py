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
from dendropy.test.support.extendedtest import ExtendedTestCase
from dendropy.test.support import datatest
from dendropy.test.support import datagen
from dendropy.test.support import pathmap
import dendropy

class TestCharStruct(ExtendedTestCase):

    def setUp(self):
        self.tb1 = dendropy.TaxonSet(label="TI1")
        for i in range(1,11):
            self.tb1.new_taxon(label="T%02d" % i)
        self.cb1 = dendropy.DnaCharacterArray(taxon_set=self.tb1, label="TI1, CA1")
        for t in self.tb1:
            self.cb1.append_taxon_sequence(t, state_symbols="AAAAAAAAAA")
        self.tb2 = dendropy.TaxonSet(label="TI2")
        for i in range(1,21):
            self.tb2.new_taxon(label="T%02d" % i)
        self.cb2 = dendropy.DnaCharacterArray(taxon_set=self.tb2, label="TI2, CA2")
        for t in self.tb2:
            self.cb2.append_taxon_sequence(t, state_symbols="CCCCCCCCCC")

class TestExtendCharacters(TestCharStruct):

    def runTest(self):
        ntax_pre = len(self.cb1)
        nchars_pre = len(self.cb1.values()[0])
        self.cb1.extend_characters(self.cb2)
        self.assertEqual(len(self.cb1), ntax_pre)
        for t in self.cb1:
            self.assertEqual(len(self.cb1[t]), 20)
            self.assertEqual(self.cb1[t].symbols_as_string(), "AAAAAAAAAACCCCCCCCCC")

class TestExtendSequencesOverwrite(TestCharStruct):

    def runTest(self):
        self.cb1.extend(self.cb2, overwrite_existing=True)
        target_ntax = 20
        self.assertEqual(len(self.cb1), target_ntax)
        self.assertEqual(len(self.cb1.taxon_set), target_ntax)
        for t in self.tb2:
            cb_tb_labels = self.cb1.taxon_set.labels()
            self.assertContained(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertContained(t.label, cb_labels)
        for t in self.cb1:
            self.assertEqual(len(self.cb1[t]), 10)
            self.assertEqual(self.cb1[t].symbols_as_string(), "CCCCCCCCCC",)

class TestExtendSequencesAppend(TestCharStruct):

    def runTest(self):
        self.cb1.extend(self.cb2, extend_existing=True)
        target_ntax = 20
        self.assertEqual(len(self.cb1), target_ntax)
        self.assertEqual(len(self.cb1.taxon_set), target_ntax)
        for t in self.tb2:
            cb_tb_labels = self.cb1.taxon_set.labels()
            self.assertContained(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertContained(t.label, cb_labels)
        for t in self.cb1:
            if int(t.label[-2:]) > 10:
                self.assertEqual(len(self.cb1[t]), 10)
                self.assertEqual(self.cb1[t].symbols_as_string(), "CCCCCCCCCC")
            else:
                self.assertEqual(len(self.cb1[t]), 20)
                self.assertEqual(self.cb1[t].symbols_as_string(), "AAAAAAAAAACCCCCCCCCC")

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
            distinct_taxa=False)

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
            distinct_taxa=False)

class CharArrayReadTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.reference_dataset = datagen.reference_single_taxonset_dataset()
        self.data_path = pathmap.named_output_path("char_read.nex")
        self.reference_dataset.write_to_path(self.data_path, "nexus")

    def testNonIndexedRead(self):
        c = dendropy.DnaCharacterArray()
        c.read_from_path(self.data_path, "nexus")
        self.assertDistinctButEqual(
            self.reference_dataset.char_arrays[0],
            c,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testIndexedRead(self):
        c = dendropy.StandardCharacterArray()
        c.read_from_stream(open(self.data_path, "rU"), "nexus", matrix_offset=1)
        self.assertDistinctButEqual(
            self.reference_dataset.char_arrays[1],
            c,
            char_type=dendropy.StandardCharacterArray,
            distinct_state_alphabets=None,
            distinct_taxa=True)

    def testIncompatibleRead(self):
        c = dendropy.DnaCharacterArray()
        self.assertRaises(ValueError, c.read_from_path, self.data_path, "nexus", matrix_offset=1)

    def testInitRead(self):
        c = dendropy.DnaCharacterArray(stream=open(self.data_path, "rU"),
                format="nexus")
        self.assertDistinctButEqual(
            self.reference_dataset.char_arrays[0],
            c,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testSameTaxaRead(self):
        c = dendropy.DnaCharacterArray()
        c.read_from_path(self.data_path,
                format="nexus",
                taxon_set=self.reference_dataset.char_arrays[0].taxon_set)
        self.assertDistinctButEqual(
            self.reference_dataset.char_arrays[0],
            c,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=False)

    def testSameTaxaInit(self):
        c = dendropy.DnaCharacterArray(stream=open(self.data_path, "rU"),
                format="nexus",
                taxon_set=self.reference_dataset.char_arrays[0].taxon_set)
        self.assertDistinctButEqual(
            self.reference_dataset.char_arrays[0],
            c,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=False)

class CharArrayWriteTest(datatest.DataObjectVerificationTestCase):

    def testDnaRountTripDistinctTaxa(self):
        c1 = datagen.reference_dna_array()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.DnaCharacterArray.get_from_path(path, "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testDnaRountTripToStringDistinctTaxa(self):
        c1 = datagen.reference_dna_array()
        s1 = c1.as_string(format="nexus")
        c2 = dendropy.DnaCharacterArray.get_from_string(s1, "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testDnaRountTripSameTaxa(self):
        c1 = datagen.reference_dna_array()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.DnaCharacterArray.get_from_path(path, "nexus", taxon_set=c1.taxon_set)
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterArray,
            distinct_state_alphabets=False,
            distinct_taxa=False)

    def testStandardRountTripDistinctTaxa(self):
        c1 = datagen.reference_standard_array()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.StandardCharacterArray.get_from_stream(open(path, "rU"), "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.StandardCharacterArray,
            distinct_state_alphabets=None,
            distinct_taxa=True)

    def testStandardRountTripSameTaxa(self):
        c1 = datagen.reference_standard_array()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.StandardCharacterArray.get_from_path(path, "nexus", taxon_set=c1.taxon_set)
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.StandardCharacterArray,
            distinct_state_alphabets=None,
            distinct_taxa=False)


if __name__ == "__main__":
    unittest.main()

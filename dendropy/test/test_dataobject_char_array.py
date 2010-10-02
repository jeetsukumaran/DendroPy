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
Tests creation, reading, update, deletion of CharMatrix objects.
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
        self.cb1 = dendropy.DnaCharacterMatrix(taxon_set=self.tb1, label="TI1, CA1")
        for t in self.tb1:
            self.cb1.append_taxon_sequence(t, state_symbols="AAAAAAAAAA")
        self.tb2 = dendropy.TaxonSet(label="TI2")
        for i in range(1,21):
            self.tb2.new_taxon(label="T%02d" % i)
        self.cb2 = dendropy.DnaCharacterMatrix(taxon_set=self.tb2, label="TI2, CA2")
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
            self.assertIn(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertIn(t.label, cb_labels)
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
            self.assertIn(t.label, cb_tb_labels)
            cb_labels = [t.label for t in self.cb1]
            self.assertIn(t.label, cb_labels)
        for t in self.cb1:
            if int(t.label[-2:]) > 10:
                self.assertEqual(len(self.cb1[t]), 10)
                self.assertEqual(self.cb1[t].symbols_as_string(), "CCCCCCCCCC")
            else:
                self.assertEqual(len(self.cb1[t]), 20)
                self.assertEqual(self.cb1[t].symbols_as_string(), "AAAAAAAAAACCCCCCCCCC")

class DnaMatrixTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.char_matrix1 = datagen.reference_dna_matrix()

    def testFromDnaCharMatrix(self):
        ca2 = dendropy.DnaCharacterMatrix(self.char_matrix1)
        self.assertDistinctButEqual(
            self.char_matrix1,
            ca2,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=False)

class StandardMatrixTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.char_matrix1 = datagen.reference_standard_matrix()

    def testFromStandardCharMatrix(self):
        ca2 = dendropy.StandardCharacterMatrix(self.char_matrix1)
        self.assertDistinctButEqual(
            self.char_matrix1,
            ca2,
            char_type=dendropy.StandardCharacterMatrix,
            distinct_state_alphabets=True,
            distinct_taxa=False)

class CharMatrixReadTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.reference_dataset = datagen.reference_single_taxonset_dataset()
        self.data_path = pathmap.named_output_path("char_read.nex")
        self.reference_dataset.write_to_path(self.data_path, "nexus")

    def testNonIndexedRead(self):
        c = dendropy.DnaCharacterMatrix()
        c.read_from_path(self.data_path, "nexus")
        self.assertDistinctButEqual(
            self.reference_dataset.char_matrices[0],
            c,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testIndexedRead(self):
        c = dendropy.StandardCharacterMatrix()
        c.read_from_stream(open(self.data_path, "rU"), "nexus", matrix_offset=1)
        self.assertDistinctButEqual(
            self.reference_dataset.char_matrices[1],
            c,
            char_type=dendropy.StandardCharacterMatrix,
            distinct_state_alphabets=None,
            distinct_taxa=True)

    def testIncompatibleRead(self):
        c = dendropy.DnaCharacterMatrix()
        self.assertRaises(ValueError, c.read_from_path, self.data_path, "nexus", matrix_offset=1)

    def testInitRead(self):
        c = dendropy.DnaCharacterMatrix(stream=open(self.data_path, "rU"),
                schema="nexus")
        self.assertDistinctButEqual(
            self.reference_dataset.char_matrices[0],
            c,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testSameTaxaRead(self):
        c = dendropy.DnaCharacterMatrix()
        c.read_from_path(self.data_path,
                schema="nexus",
                taxon_set=self.reference_dataset.char_matrices[0].taxon_set)
        self.assertDistinctButEqual(
            self.reference_dataset.char_matrices[0],
            c,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=False)

    def testSameTaxaInit(self):
        c = dendropy.DnaCharacterMatrix(stream=open(self.data_path, "rU"),
                schema="nexus",
                taxon_set=self.reference_dataset.char_matrices[0].taxon_set)
        self.assertDistinctButEqual(
            self.reference_dataset.char_matrices[0],
            c,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=False)

class CharMatrixWriteTest(datatest.DataObjectVerificationTestCase):

    def testDnaRountTripDistinctTaxa(self):
        c1 = datagen.reference_dna_matrix()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.DnaCharacterMatrix.get_from_path(path, "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testDnaRountTripToStringDistinctTaxa(self):
        c1 = datagen.reference_dna_matrix()
        s1 = c1.as_string(schema="nexus")
        c2 = dendropy.DnaCharacterMatrix.get_from_string(s1, "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=True)

    def testDnaRountTripSameTaxa(self):
        c1 = datagen.reference_dna_matrix()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.DnaCharacterMatrix.get_from_path(path, "nexus", taxon_set=c1.taxon_set)
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.DnaCharacterMatrix,
            distinct_state_alphabets=False,
            distinct_taxa=False)

    def testStandardRountTripDistinctTaxa(self):
        c1 = datagen.reference_standard_matrix()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.StandardCharacterMatrix.get_from_stream(open(path, "rU"), "nexus")
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.StandardCharacterMatrix,
            distinct_state_alphabets=None,
            distinct_taxa=True)

    def testStandardRountTripSameTaxa(self):
        c1 = datagen.reference_standard_matrix()
        path = pathmap.named_output_path("char_rw_dna.nex")
        c1.write_to_path(path, "nexus")
        c2 = dendropy.StandardCharacterMatrix.get_from_path(path, "nexus", taxon_set=c1.taxon_set)
        self.assertDistinctButEqual(
            c1,
            c2,
            char_type=dendropy.StandardCharacterMatrix,
            distinct_state_alphabets=None,
            distinct_taxa=False)


if __name__ == "__main__":
    unittest.main()

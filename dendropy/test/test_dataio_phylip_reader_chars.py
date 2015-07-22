# !/usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests for general PHYLIP character matrix reading.
"""

import unittest
import dendropy
import collections
from dendropy.utility import error
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap
from dendropy.test.support import standard_file_test_chars
from dendropy.test.support import compare_and_validate
from dendropy.dataio import phylipreader
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class PhylipCharactersReaderDnaTestCase(
        standard_file_test_chars.DnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_phylip(self):
        src_filenames = [
                "standard-test-chars-dna.relaxed.phylip",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.DnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="phylip",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class PhylipCharactersReaderRnaTestCase(
        standard_file_test_chars.RnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_phylip(self):
        src_filenames = [
                "standard-test-chars-rna.relaxed.phylip",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.RnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="phylip",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class PhylipCharactersReaderProteinTestCase(
        standard_file_test_chars.ProteinTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_phylip(self):
        src_filenames = [
                "standard-test-chars-protein.relaxed.phylip",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ProteinCharacterMatrix,
                    src_filepath=src_path,
                    schema="phylip",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class PhylipStandardCharacters01234TestCase(
        standard_file_test_chars.Standard01234TestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build(state_alphabet_fundamental_symbols="0123456789-")

    def test_basic_phylip(self):
        src_filenames = [
                "standard-test-chars-generic.relaxed.phylip",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.StandardCharacterMatrix,
                    src_filepath=src_path,
                    schema="phylip",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class PhylipCharactersContinuousTestCase(
        standard_file_test_chars.ContinuousTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                ("standard-test-chars-continuous.relaxed.phylip", {}),
                ("standard-test-chars-continuous.interleaved.phylip", {"interleaved": True}),
                ]
        for src_idx, (src_filename, kwargs) in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ContinuousCharacterMatrix,
                    src_filepath=src_path,
                    schema="phylip",
                    factory_kwargs=kwargs,
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class PhylipVariantsTestCases(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.expected_seqs = collections.OrderedDict()
        cls.expected_seqs["Turkey"]     = "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT"
        cls.expected_seqs["Salmo gair"] = "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT"
        cls.expected_seqs["H. Sapiens"] = "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA"
        cls.expected_seqs["Chimp"]      = "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT"
        cls.expected_seqs["Gorilla"]    = "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA"

    def test_strict_sequential(self):
        s = """\
  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT
GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT
ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT
AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT
AAACCATTGC CGGTACGCTT AA
        """
        char_matrix = dendropy.DnaCharacterMatrix.get_from_string(
                s,
                "phylip",
                strict=True)
        self.assertEqual(len(char_matrix), len(self.expected_seqs))
        self.assertEqual(len(char_matrix.taxon_namespace), len(self.expected_seqs))
        for taxon, expected_taxon in zip(char_matrix, self.expected_seqs):
            self.assertEqual(taxon.label, expected_taxon)
            self.assertEqual(char_matrix[taxon].symbols_as_string(), self.expected_seqs[expected_taxon])

    def test_strict_interleaved(self):
        s = """\
5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
Salmo gairAAGCCTTGGC AGTGCAGGGT
H. SapiensACCGGTTGGC CGTTCAGGGT
Chimp     AAACCCTTGC CGTTACGCTT
Gorilla   AAACCCTTGC CGGTACGCTT

GAGCCCGGGC AATACAGGGT AT
GAGCCGTGGC CGGGCACGGT AT
ACAGGTTGGC CGTTCAGGGT AA
AAACCGAGGC CGGGACACTC AT
AAACCATTGC CGGTACGCTT AA
        """
        char_matrix = dendropy.DnaCharacterMatrix.get_from_string(
                s,
                "phylip",
                interleaved=True,
                strict=True)
        self.assertEqual(len(char_matrix), len(self.expected_seqs))
        self.assertEqual(len(char_matrix.taxon_namespace), len(self.expected_seqs))
        for taxon, expected_taxon in zip(char_matrix, self.expected_seqs):
            self.assertEqual(taxon.label, expected_taxon)
            self.assertEqual(char_matrix[taxon].symbols_as_string(), self.expected_seqs[expected_taxon])

    def test_strict_interleaved_with_bad_chars(self):
        s = """\
5    42
Turkey    AAGCTNGGGC ATTTCA3828GGGT
Salmo gairAAGCCTTGGC AGTGCA3828GGGT
H. SapiensACCGGTTGGC CGTTCA3828GGGT
Chimp     AAACCCTTGC CGTTAC3828GCTT
Gorilla   AAACCCTTGC CGGTAC3828GCTT

GAGCCCGGGC AATACAGGGT AT
GAGCCGTGGC CGGGCACGGT AT
ACAGGTTGGC CGTTCAGGGT AA
AAACCGAGGC CGGGACACTC AT
AAACCATTGC CGGTACGCTT AA
        """
        char_matrix = dendropy.DnaCharacterMatrix.get_from_string(
                s,
                "phylip",
                interleaved=True,
                strict=True,
                ignore_invalid_chars=True)
        self.assertEqual(len(char_matrix), len(self.expected_seqs))
        self.assertEqual(len(char_matrix.taxon_namespace), len(self.expected_seqs))
        for taxon, expected_taxon in zip(char_matrix, self.expected_seqs):
            self.assertEqual(taxon.label, expected_taxon)
            self.assertEqual(char_matrix[taxon].symbols_as_string(), self.expected_seqs[expected_taxon])

class PhylipContinuousVariantsTestCases(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.expected_seqs = collections.OrderedDict()
        cls.expected_seqs["Turkey"]     = [-231.6391 ,  972.4189  ,  626.6717  ,  -328.6811 ,  -213.5738 ,  464.3897  ,  -91.3483  ,  349.8176  ,  333.4800  ,  521.4970  ]
        cls.expected_seqs["Salmo gair"] = [ 104.4199 ,  669.7402  ,  -68.6082  ,  975.4302  ,  -874.4510 ,  -191.3305 ,  -179.8437 ,  655.5611  ,  -657.4532 ,  -563.7863 ]
        cls.expected_seqs["H. Sapiens"] = [-613.2947 ,  -600.7053 ,  -700.5140 ,  438.6092  ,  615.5268  ,  640.7933  ,  503.8948  ,  -159.7922 ,  866.8036  ,  274.0275  ]
        cls.expected_seqs["Chimp"]      = [-654.7695 ,  103.3806  ,  -971.8866 ,  853.9164  ,  653.5797  ,  823.6672  ,  -476.6859 ,  325.9331  ,  456.0902  ,  -399.7095 ]
        cls.expected_seqs["Gorilla"]    = [-762.4904 ,  808.3665  ,  522.5775  ,  250.6523  ,  -287.9786 ,  -995.4612 ,  571.9263  ,  -793.3975 ,  -42.7027  ,  186.8869  ]

    def test_strict_sequential(self):
        s = """\
  5    10
Turkey    -231.6391 972.4189  626.6717  -328.6811 -213.5738 464.3897  -91.3483  349.8176  333.4800  521.4970
Salmo gair 104.4199 669.7402  -68.6082  975.4302  -874.4510 -191.3305 -179.8437 655.5611  -657.4532 -563.7863
H. Sapiens-613.2947 -600.7053 -700.5140 438.6092  615.5268  640.7933  503.8948  -159.7922 866.8036  274.0275
Chimp     -654.7695 103.3806  -971.8866 853.9164  653.5797  823.6672  -476.6859 325.9331  456.0902  -399.7095
Gorilla   -762.4904 808.3665  522.5775  250.6523  -287.9786 -995.4612 571.9263  -793.3975 -42.7027  186.8869
        """
        char_matrix = dendropy.ContinuousCharacterMatrix.get_from_string(
                s,
                "phylip",
                strict=True)
        self.assertEqual(len(char_matrix), len(self.expected_seqs))
        self.assertEqual(len(char_matrix.taxon_namespace), len(self.expected_seqs))
        for taxon, expected_taxon in zip(char_matrix, self.expected_seqs):
            self.assertEqual(taxon.label, expected_taxon)
            self.assertEqual(char_matrix[taxon].values(), self.expected_seqs[expected_taxon])

    def test_strict_interleaved(self):
        s = """\
  5    10
Turkey    -231.6391 972.4189  626.6717  -328.6811
Salmo gair 104.4199 669.7402  -68.6082  975.4302
H. Sapiens-613.2947 -600.7053 -700.5140 438.6092
Chimp     -654.7695 103.3806  -971.8866 853.9164
Gorilla   -762.4904 808.3665  522.5775  250.6523
 -213.5738 464.3897  -91.3483  349.8176  333.4800  521.4970
 -874.4510 -191.3305 -179.8437 655.5611  -657.4532 -563.7863
 615.5268  640.7933  503.8948  -159.7922 866.8036  274.0275
 653.5797  823.6672  -476.6859 325.9331  456.0902  -399.7095
 -287.9786 -995.4612 571.9263  -793.3975 -42.7027  186.8869
        """
        char_matrix = dendropy.ContinuousCharacterMatrix.get_from_string(
                s,
                "phylip",
                interleaved=True,
                strict=True)
        self.assertEqual(len(char_matrix), len(self.expected_seqs))
        self.assertEqual(len(char_matrix.taxon_namespace), len(self.expected_seqs))
        for taxon, expected_taxon in zip(char_matrix, self.expected_seqs):
            self.assertEqual(taxon.label, expected_taxon)
            self.assertEqual(char_matrix[taxon].values(), self.expected_seqs[expected_taxon])

if __name__ == "__main__":
    unittest.main()

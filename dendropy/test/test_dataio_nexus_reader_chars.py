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
Tests for general NEXUS character matrix reading.
"""

import unittest
import dendropy
from dendropy.utility import error
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap
from dendropy.test.support import standard_file_test_chars
from dendropy.test.support import compare_and_validate
from dendropy.dataio import nexusreader
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class NexusCharactersReaderDnaTestCase(
        standard_file_test_chars.DnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                "standard-test-chars-dna.simple.nexus",
                "standard-test-chars-dna.basic.nexus",
                "standard-test-chars-dna.interleaved.nexus",
                "standard-test-chars-dna.matchchar.nexus",
                "standard-test-chars-dna.multi.nexus",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.DnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexus",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexusCharactersReaderRnaTestCase(
        standard_file_test_chars.RnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                "standard-test-chars-rna.simple.nexus",
                "standard-test-chars-rna.basic.nexus",
                "standard-test-chars-rna.interleaved.nexus",
                "standard-test-chars-rna.matchchar.nexus",
                "standard-test-chars-rna.multi.nexus",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.RnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexus",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexusCharactersReaderProteinTestCase(
        standard_file_test_chars.ProteinTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                "standard-test-chars-protein.simple.nexus",
                "standard-test-chars-protein.basic.nexus",
                "standard-test-chars-protein.interleaved.nexus",
                "standard-test-chars-protein.matchchar.nexus",
                "standard-test-chars-protein.multi.nexus",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ProteinCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexus",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexusCharactersContinuousTestCase(
        standard_file_test_chars.ContinuousTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                "standard-test-chars-continuous.mesquite.nexus",
                "standard-test-chars-continuous.mesquite.interleaved.nexus",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ContinuousCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexus",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexusStandardCharacters01234TestCase(
        standard_file_test_chars.Standard01234TestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexus(self):
        src_filenames = [
                "standard-test-chars-generic.simple.nexus",
                "standard-test-chars-generic.basic.nexus",
                "standard-test-chars-generic.dotted.nexus",
                "standard-test-chars-generic.interleaved.nexus",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.StandardCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexus",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)


class NexusTooManyTaxaTest(
        dendropytest.ExtendedTestCase):

    def testTooManyTaxaNonInterleaved(self):
        data_str = """\
        #NEXUS
        BEGIN TAXA;
            DIMENSIONS NTAX=2;
            TAXLABELS AAA BBB ;
        END;
        BEGIN CHARACTERS;
            DIMENSIONS  NCHAR=8;
            FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
            MATRIX
                AAA ACGTACGT
                BBB ACGTACGT
                CCC ACGTACGT
            ;
        END;
        """
        self.assertRaises(nexusreader.NexusReader.TooManyTaxaError,
                dendropy.DnaCharacterMatrix.get_from_string,
                data_str,
                'nexus')

class NexusCharsSubsetsTest(
        compare_and_validate.Comparator,
        dendropytest.ExtendedTestCase):

    def verify_subsets(self, src_filename, expected_sets):
        """
        ``src_filename`` -- name of file containing full data and charsets
                          statement
        ``expected_sets`` -- dictionary with keys = label of charset, and values
                           = name of file with subset of characters correspond
                           to the charset.
        """

        src_data = dendropy.DnaCharacterMatrix.get_from_path(
                pathmap.char_source_path(src_filename),
                'nexus')

        state_alphabet = src_data.default_state_alphabet
        self.assertEqual(len(src_data.character_subsets), len(expected_sets))
        for label, expected_data_file in expected_sets.items():

            _LOG.debug(label)

            self.assertTrue(label in src_data.character_subsets)
            result_subset = src_data.export_character_subset(label)
            expected_subset = dendropy.DnaCharacterMatrix.get_from_path(
                pathmap.char_source_path(expected_data_file),
                'nexus')

            # confirm subset is correct
            self.compare_distinct_char_matrix(
                    result_subset,
                    expected_subset,
                    taxon_namespace_scoped=False,
                    )

            # mutate new and confirm that old remains unchanged
            e1_symbols = src_data[0].symbols_as_string()
            r1 = result_subset[0]
            dummy_state = state_alphabet["A"]
            for idx in range(len(r1)):
                r1[idx].value = dummy_state
            self.assertEqual(e1_symbols, src_data[0].symbols_as_string())

            # mutate old and confirm that new remains unchanged
            r2_symbols = result_subset[1].symbols_as_string()
            e2 = src_data[1]
            dummy_state = state_alphabet["A"]
            for idx in range(len(e2)):
                e2[idx].value = dummy_state
            self.assertEqual(r2_symbols, result_subset[1].symbols_as_string())

    def testNonInterleaved(self):
        """
        Charsets here go through all forms of position specification.
        """
        expected_sets = {
                "coding" : "primates.chars.subsets-coding.nexus",
                "noncoding" : "primates.chars.subsets-noncoding.nexus",
                "1stpos" : "primates.chars.subsets-1stpos.nexus",
                "2ndpos" : "primates.chars.subsets-2ndpos.nexus",
                "3rdpos" : "primates.chars.subsets-3rdpos.nexus",
                }
        self.verify_subsets('primates.chars.subsets-all.nexus', expected_sets)

    def testInterleaved(self):
        """
        A bug in DendroPy resulted in the block immediately following an
        interleaved character matrix DATA or CHARACTERS block being skipped.
        This tests for it by ensuring that the ASSUMPTIONS block following an
        interleaved CHARACTERS block is parsed. A better test would approach
        the issue more directly, by checking to see if block parsing left the
        stream reader in the correct position.
        """
        expected_sets = {
                "c1" : "interleaved-charsets-c1.nex",
                "c2" : "interleaved-charsets-c2.nex",
                "c3" : "interleaved-charsets-c3.nex",
                }
        self.verify_subsets('interleaved-charsets-all.nex', expected_sets)

if __name__ == "__main__":
    unittest.main()

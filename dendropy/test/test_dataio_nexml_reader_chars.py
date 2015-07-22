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
from dendropy.dataio import nexmlreader
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class NexmlCharactersReaderDnaTestCase(
        standard_file_test_chars.DnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexml(self):
        src_filenames = [
                "standard-test-chars-dna.as_cells.nexml",
                "standard-test-chars-dna.as_seqs.nexml",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.DnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexml",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexmlCharactersReaderRnaTestCase(
        standard_file_test_chars.RnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexml(self):
        src_filenames = [
                "standard-test-chars-rna.as_cells.nexml",
                "standard-test-chars-rna.as_seqs.nexml",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.RnaCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexml",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexmlCharactersReaderProteinTestCase(
        standard_file_test_chars.ProteinTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexml(self):
        src_filenames = [
                "standard-test-chars-protein.as_cells.nexml",
                "standard-test-chars-protein.as_seqs.nexml",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ProteinCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexml",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexmlCharactersContinuousTestCase(
        standard_file_test_chars.ContinuousTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexml(self):
        src_filenames = [
                "standard-test-chars-continuous.as_cells.nexml",
                "standard-test-chars-continuous.as_seqs.nexml",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.ContinuousCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexml",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

class NexmlStandardCharacters01234TestCase(
        standard_file_test_chars.Standard01234TestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_nexml(self):
        src_filenames = [
                "standard-test-chars-generic.as_cells.nexml",
                "standard-test-chars-generic.as_seqs.nexml",
                ]
        for src_idx, src_filename in enumerate(src_filenames):
            # print(src_idx, src_filename)
            src_path = pathmap.char_source_path(src_filename)
            self.verify_get_from(
                    matrix_type=dendropy.StandardCharacterMatrix,
                    src_filepath=src_path,
                    schema="nexml",
                    factory_kwargs={},
                    check_taxon_annotations=False,
                    check_matrix_annotations=False,
                    check_sequence_annotations=False,
                    check_column_annotations=False,
                    check_cell_annotations=False)

if __name__ == "__main__":
    unittest.main()

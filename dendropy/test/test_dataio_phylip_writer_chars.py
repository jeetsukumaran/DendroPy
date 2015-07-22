#! /usr/bin/env python

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
Tests for PHYLIP tree list writing.
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import compare_and_validate
from dendropy.test.support import pathmap
from dendropy.test.support import standard_file_test_chars

class PhylipWriterCharactersTestCase(
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.check_taxon_annotations = False
        cls.check_matrix_annotations = False
        cls.check_sequence_annotations = False
        cls.check_column_annotations = False
        cls.check_cell_annotations = False
        standard_file_test_chars.DnaTestChecker.build()
        standard_file_test_chars.RnaTestChecker.build()
        standard_file_test_chars.ProteinTestChecker.build()
        standard_file_test_chars.Standard01234TestChecker.build()
        standard_file_test_chars.ContinuousTestChecker.build()
        cls.srcs = (
                ("standard-test-chars-dna.relaxed.phylip", dendropy.DnaCharacterMatrix, standard_file_test_chars.DnaTestChecker),
                ("standard-test-chars-rna.relaxed.phylip", dendropy.RnaCharacterMatrix, standard_file_test_chars.RnaTestChecker),
                ("standard-test-chars-protein.relaxed.phylip", dendropy.ProteinCharacterMatrix, standard_file_test_chars.ProteinTestChecker),
                ("standard-test-chars-generic.relaxed.phylip", dendropy.StandardCharacterMatrix, standard_file_test_chars.Standard01234TestChecker),
                ("standard-test-chars-continuous.relaxed.phylip", dendropy.ContinuousCharacterMatrix, standard_file_test_chars.ContinuousTestChecker),
                )

    def verify_char_matrix(self, char_matrix, src_matrix_checker_type):
        self.assertEqual(type(char_matrix), src_matrix_checker_type.matrix_type)
        if src_matrix_checker_type.matrix_type is dendropy.StandardCharacterMatrix:
            src_matrix_checker_type.create_class_fixtures_label_sequence_map_based_on_state_alphabet(src_matrix_checker_type,
                    char_matrix.default_state_alphabet)
        standard_file_test_chars.general_char_matrix_checker(
                self,
                char_matrix,
                src_matrix_checker_type,
                check_taxon_annotations=self.check_taxon_annotations,
                check_matrix_annotations=self.check_matrix_annotations,
                check_sequence_annotations=self.check_sequence_annotations,
                check_column_annotations=self.check_column_annotations,
                check_cell_annotations=self.check_cell_annotations,)

    def test_basic_phylip_chars(self):
        for src_filename, matrix_type, src_matrix_checker_type in self.__class__.srcs:
            src_path = pathmap.char_source_path(src_filename)
            d1 = matrix_type.get_from_path(src_path, "phylip")
            for strict in (True, False):
                for spaces_to_underscores in (True, False):
                    for force_unique_taxon_labels in (True, False):
                        s = self.write_out_validate_equal_and_return(
                                d1, "phylip", {
                                    "strict": strict,
                                    "spaces_to_underscores" : spaces_to_underscores,
                                    "force_unique_taxon_labels" : force_unique_taxon_labels,
                                    })
                        d2 = matrix_type.get_from_string(s, "phylip")
                        self.verify_char_matrix(d2, src_matrix_checker_type)

class PhylipWriterCharactersVariantsTestCase(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        d = dendropy.DnaCharacterMatrix()
        cls.original_labels = [
                "a0_123456789_1",
                "a0_123456789_2",
                "a0_123456789_3",
                "a0_123456789_4",
                "a0_123456789_5",
                "b0_123456789_1",
                "b0_123456789_2",
                "b0_123456789_3",
                "b0_123456789_4",
                ]
        for label in cls.original_labels:
            t = d.taxon_namespace.require_taxon(label=label)
            d[t] = d.default_state_alphabet.get_states_for_symbols("AACGT")
        cls.data = d

    def test_strict_write(self):
        s0 = self.data.as_string("phylip", strict=True)
        d2 = dendropy.DnaCharacterMatrix.get_from_string(s0, "phylip", strict=True)
        obs_labels = set([t.label for t in d2])
        self.assertEqual(len(obs_labels), len(self.original_labels))
        for t in d2:
            self.assertEqual(len(t.label), 10)
            self.assertEqual(d2[t].symbols_as_string(), "AACGT")

if __name__ == "__main__":
    unittest.main()

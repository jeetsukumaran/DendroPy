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
Tests for general FASTA reading.
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap
from dendropy.test.support import standard_file_test_chars

class FastaDnaReaderTestCase(
        standard_file_test_chars.DnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_fasta(self):
        src_path = pathmap.char_source_path("standard-test-chars-dna.fasta")
        self.verify_get_from(
                matrix_type=dendropy.DnaCharacterMatrix,
                src_filepath=src_path,
                schema="fasta",
                factory_kwargs={},
                check_taxon_annotations=False,
                check_matrix_annotations=False,
                check_sequence_annotations=False,
                check_column_annotations=False,
                check_cell_annotations=False)

class FastaRnaReaderTestCase(
        standard_file_test_chars.RnaTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_fasta(self):
        src_path = pathmap.char_source_path("standard-test-chars-rna.fasta")
        self.verify_get_from(
                matrix_type=dendropy.RnaCharacterMatrix,
                src_filepath=src_path,
                schema="fasta",
                factory_kwargs={},
                check_taxon_annotations=False,
                check_matrix_annotations=False,
                check_sequence_annotations=False,
                check_column_annotations=False,
                check_cell_annotations=False)

class FastaProteinReaderTestCase(
        standard_file_test_chars.ProteinTestChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.build()

    def test_basic_fasta(self):
        src_path = pathmap.char_source_path("standard-test-chars-protein.fasta")
        self.verify_get_from(
                matrix_type=dendropy.ProteinCharacterMatrix,
                src_filepath=src_path,
                schema="fasta",
                factory_kwargs={},
                check_taxon_annotations=False,
                check_matrix_annotations=False,
                check_sequence_annotations=False,
                check_column_annotations=False,
                check_cell_annotations=False)

if __name__ == "__main__":
    unittest.main()

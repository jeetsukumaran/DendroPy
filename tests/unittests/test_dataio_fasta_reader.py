# !/usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Tests for general FASTA reading.
"""

import unittest
import dendropy
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest
from support import pathmap
from support import standard_file_test_chars

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

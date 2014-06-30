#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Tests for NEXUS tree list writing.
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import compare_and_validate
from dendropy.test.support import pathmap
from dendropy.test.support import standard_file_test_chars

class NexusWriterCharactersTestCase(
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
        cls.srcs = (
                ("standard-test-chars-dna.multi.nexus", standard_file_test_chars.DnaTestChecker),
                ("standard-test-chars-rna.multi.nexus", standard_file_test_chars.RnaTestChecker),
                ("standard-test-chars-protein.multi.nexus", standard_file_test_chars.ProteinTestChecker),
                ("standard-test-chars-generic.interleaved.nexus", standard_file_test_chars.Standard01234TestChecker),
                )

    def test_basic_nexus_chars_dna(self):
        src_filename = "standard-test-chars-dna.basic.nexus"
        src_path = pathmap.char_source_path(src_filename)
        d1 = dendropy.DnaCharacterMatrix.get_from_path(
                src_path,
                "nexus")
        s = self.write_out_validate_equal_and_return(
                d1, "nexus", {})
        d2 = dendropy.DnaCharacterMatrix.get_from_string(
                s,
                "nexus")

if __name__ == "__main__":
    unittest.main()

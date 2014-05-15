# !/usr/bin/env python

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
Tests for general NEXUS character matrix reading.
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap
from dendropy.test.support import datagen_standard_file_test_chars

class NexusCharactersReaderTestCase(
        datagen_standard_file_test_chars.DnaTestChecker,
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

if __name__ == "__main__":
    unittest.main()

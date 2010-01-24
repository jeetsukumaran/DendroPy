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
Phylip data read/write parse/format tests.
"""

import sys
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import extendedtest
from dendropy.dataio import phylip

strict_dna_valid_chars_clean_str = """\
  5    20
Turkey    AAGCTNGGGC ATTTCAGGGT
Salmo gairAAGCCTTGGC AGTGCAGGGT
H. SapiensACCGGTTGGC CGTTCAGGGT
Chimp     AAACCCTTGC CGTTACGCTT
Gorilla   AAACCCTTGC CGGTACGCTT
"""

strict_dna_valid_chars_messy_str = "\n".join([
    "  5    20    ",
    "Turkey    AAGCTNGGGC ATTTCAGGGT    ",
    "Salmo gair     ",
    "  ",
    "",
    "    AAGCCTTG   GC  ",
    "AGTGCAG   ",
    "GGT",
    "H. SapiensA",
    "CCGGTTGGC CGTTCAGGGT",
    "Chimp     ",
    "AAACCCTTGC ",
    "   CGTTACGCTT",
    "Gorilla   AAACCCTTGC CGGTACGCTT",
    "          ",
    "          ",
    ])

strict_dna_invalid_chars_messy_str = "\n".join([
    "  5    20    ",
    "Turkey    AAGCT441N   z q lGGGC ATTTCAGGGT    ",
    "Salmo gair     ",
    "  ",
    "",
    "    AAGCCTTG   GC  ",
    "AGTGCAG   ",
    "GGT",
    "H. SapiensA",
    "CCGGTTGGC CGTTCAGGq   *. q.u qGT",
    "Chimp     ",
    "AAACCCTT%56%qeq qGC ",
    "   CGTTACGCTT",
    "Gorilla   AAACC33CTTGC CGGTACGCTT",
    "          ",
    "          ",
    ])

strict_dna_expected_taxon_labels = ["Turkey", "Salmo gair", "H. Sapiens", "Chimp", "Gorilla"]
strict_dna_expected_seq_symbols = [
        "AAGCTNGGGCATTTCAGGGT",
        "AAGCCTTGGCAGTGCAGGGT",
        "ACCGGTTGGCCGTTCAGGGT",
        "AAACCCTTGCCGTTACGCTT",
        "AAACCCTTGCCGGTACGCTT",
]

class StrictSequentialDnaTest(extendedtest.ExtendedTestCase):

    def verify(self, dataset):
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqual(len(dataset.tree_lists), 0)
        taxon_set = dataset.taxon_sets[0]
        char_matrix = dataset.char_matrices[0]
        self.assertIs(taxon_set, char_matrix.taxon_set)
        self.assertEqual(len(char_matrix), len(strict_dna_expected_seq_symbols))
        for i, taxon in enumerate(taxon_set):
            self.assertIn(taxon, char_matrix)
            self.assertIs(char_matrix[taxon], char_matrix[i])
            seqs = char_matrix[taxon].symbols_as_string()
            self.assertEqual(seqs, strict_dna_expected_seq_symbols[i])

    def test_strict_dna_valid_chars_clean_str(self):
        pr = phylip.PhylipReader(strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        dataset = pr.read(StringIO(strict_dna_valid_chars_clean_str))
        self.verify(dataset)

    def test_strict_dna_valid_chars_messy_str(self):
        pr = phylip.PhylipReader(strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        dataset = pr.read(StringIO(strict_dna_valid_chars_messy_str))
        self.verify(dataset)

    def test_strict_dna_invalid_chars_messy_str(self):
        pr = phylip.PhylipReader(strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=True)
        dataset = pr.read(StringIO(strict_dna_invalid_chars_messy_str))
        self.verify(dataset)

    def test_raises_error_strict_dna_invalid_chars_messy_str(self):
        pr = phylip.PhylipReader(strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        self.assertRaises(phylip.PhylipReader.PhylipStrictSequentialError,
            pr.read,
            StringIO(strict_dna_invalid_chars_messy_str))

if __name__ == "__main__":
    unittest.main()

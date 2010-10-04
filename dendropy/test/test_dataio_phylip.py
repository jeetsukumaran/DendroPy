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
Phylip data read/write parse/format tests.
"""

import sys
import unittest
import tempfile
from cStringIO import StringIO

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import extendedtest
from dendropy.test.support import datatest
from dendropy.dataio import phylip
import dendropy

###############################################################################
## STRICT SEQUENTIAL DATA

strict_sequential_dna_valid_chars_clean_str = """\
  5    20
Turkey    AAGCTNGGGC ATTTCAGGGT
Salmo gairAAGCCTTGGC AGTGCAGGGT
H. SapiensACCGGTTGGC CGTTCAGGGT
Chimp     AAACCCTTGC CGTTACGCTT
Gorilla   AAACCCTTGC CGGTACGCTT
"""

strict_sequential_dna_valid_chars_messy_str = "\n".join([
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

strict_sequential_dna_invalid_chars_messy_str = "\n".join([
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

strict_sequential_dna_expected_taxon_labels = [
        "Turkey",
        "Salmo gair",
        "H. Sapiens",
        "Chimp",
        "Gorilla"
]

strict_sequential_dna_expected_seq_symbols = [
        "AAGCTNGGGCATTTCAGGGT",
        "AAGCCTTGGCAGTGCAGGGT",
        "ACCGGTTGGCCGTTCAGGGT",
        "AAACCCTTGCCGTTACGCTT",
        "AAACCCTTGCCGGTACGCTT",
]

###############################################################################
## RELAXED SEQUENTIAL DATA

relaxed_sequential_dna_valid_chars_clean_str = """\
  5    20
Taxon1_ABCDEFG$AAGCTNGGGC ATTTCAGGGT
Taxon2_ABCDEFG$AAGCCTTGGC AGTGCAGGGT
Taxon3_ABCDEFG$ACCGGTTGGC CGTTCAGGGT
Taxon4_ABCDEFG$AAACCCTTGC CGTTACGCTT
Taxon5_ABCDEFG$AAACCCTTGC CGGTACGCTT
"""

relaxed_sequential_dna_valid_chars_messy_str = "\n".join([
    "  5    20    ",
    "Taxon1_ABCDEFG$AAGCTNGGGC ATTTCAGGGT    ",
    "Taxon2_ABCDEFG     ",
    "  ",
    "",
    "    AAGCCTTG   GC  ",
    "AGTGCAG   ",
    "GGT",
    "Taxon3_ABCDEFG$A",
    "CCGGTTGGC CGTTCAGGGT",
    "Taxon4_ABCDEFG     ",
    "AAACCCTTGC ",
    "   CGTTACGCTT",
    "Taxon5_ABCDEFG$AAACCCTTGC CGGTACGCTT",
    "          ",
    "          ",
    ])

relaxed_sequential_dna_invalid_chars_messy_str = "\n".join([
    "  5    20    ",
    "Taxon1_ABCDEFG$AAGCTNq q q z 44 212 1GGGC ATTTCA  123 1GGGT    ",
    "Taxon2_ABCDEFG$112     ",
    "  41 121 1",
    "",
    "    AAGCC412341 14 141 1TTG    1 GC  ",
    "AGTG4 411 1 1 1 CAG   ",
    "GGT4",
    "Taxon3_ABCDEFG$A",
    "CCG414 141 GTTGGC CGTTCAGGGT",
    "Taxon4_ABCDEFG     331",
    "AAACCCT333TGC ",
    "   CGTTAq3 q3 34CGCTT",
    "Taxon5_ABCDEFG$AAA 3431 CCCTTGC CGGTACGCTT",
    "           ",
    "           ",
    ])

relaxed_sequential_dna_expected_taxon_labels = [
        "Taxon1_ABCDEFG",
        "Taxon2_ABCDEFG",
        "Taxon3_ABCDEFG",
        "Taxon4_ABCDEFG",
        "Taxon5_ABCDEFG"
]

relaxed_sequential_dna_expected_seq_symbols = [
        "AAGCTNGGGCATTTCAGGGT",
        "AAGCCTTGGCAGTGCAGGGT",
        "ACCGGTTGGCCGTTCAGGGT",
        "AAACCCTTGCCGTTACGCTT",
        "AAACCCTTGCCGGTACGCTT",
]

###############################################################################
## STRICT INTERLEAVED DATA

strict_interleaved_dna_valid_chars_clean_str = """\
    5 100

Lagomorpha?????????? ?????????C CAATCTACAC ACGGG-GTAG GGATTACATA
Human     AGCCACACCC TAGGGTTGGC CAATCTACTC CCAGGAGCAG GGAGGGCAGG
Opossum   AGCCACACCC CAACCTTAGC CAATAGACAT CCAGAAGCCC AAAAGGCAAG
Chicken   GCCCGGGGAA GAGGAGGGGC CCGGCGG-AG GCGATAAAAG TGGGGACACA
Frog      GGATGGAGAA TTAGAGCACT TGTTCTTTTT GCAGAAGCTC AGAATAAACG

TTTGGATGGT AG---GATAT GGGCCTACCA TGGCGTTAAC GGGT-AACGY
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGY
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
TTTCGATGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
"""

strict_interleaved_dna_expected_taxon_labels = [
        "Lagomorpha",
        "Human",
        "Opossum",
        "Chicken",
        "Frog"
]

strict_interleaved_dna_expected_seq_symbols = [
        "???????????????????CCAATCTACACACGGG-GTAGGGATTACATATTTGGATGGTAG---GATATGGGCCTACCATGGCGTTAACGGGT-AACGY",
        "AGCCACACCCTAGGGTTGGCCAATCTACTCCCAGGAGCAGGGAGGGCAGGTTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
        "AGCCACACCCCAACCTTAGCCAATAGACATCCAGAAGCCCAAAAGGCAAGTTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGY",
        "GCCCGGGGAAGAGGAGGGGCCCGGCGG-AGGCGATAAAAGTGGGGACACATTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
        "GGATGGAGAATTAGAGCACTTGTTCTTTTTGCAGAAGCTCAGAATAAACGTTTCGATGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
]

###############################################################################
## STRICT INTERLEAVED DATA

relaxed_interleaved_dna_valid_chars_clean_str = """\
    5 100

Taxon1_ABCDEFG$?????????? ?????????C CAATCTACAC ACGGG-GTAG GGATTACATA
Taxon2_ABCDEFG$AGCCACACCC TAGGGTTGGC CAATCTACTC CCAGGAGCAG GGAGGGCAGG
Taxon3_ABCDEFG$AGCCACACCC CAACCTTAGC CAATAGACAT CCAGAAGCCC AAAAGGCAAG
Taxon4_ABCDEFG$GCCCGGGGAA GAGGAGGGGC CCGGCGG-AG GCGATAAAAG TGGGGACACA
Taxon5_ABCDEFG$GGATGGAGAA TTAGAGCACT TGTTCTTTTT GCAGAAGCTC AGAATAAACG

TTTGGATGGT AG---GATAT GGGCCTACCA TGGCGTTAAC GGGT-AACGY
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGY
TTTCGACGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
TTTCGATGGT AA---GGTAT TGGCTTACCG TGGCAATGAC AGGT-GACGG
"""

relaxed_interleaved_dna_expected_taxon_labels = [
        "Taxon1_ABCDEFG",
        "Taxon2_ABCDEFG",
        "Taxon3_ABCDEFG",
        "Taxon4_ABCDEFG",
        "Taxon5_ABCDEFG"
]

relaxed_interleaved_dna_expected_seq_symbols = [
        "???????????????????CCAATCTACACACGGG-GTAGGGATTACATATTTGGATGGTAG---GATATGGGCCTACCATGGCGTTAACGGGT-AACGY",
        "AGCCACACCCTAGGGTTGGCCAATCTACTCCCAGGAGCAGGGAGGGCAGGTTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
        "AGCCACACCCCAACCTTAGCCAATAGACATCCAGAAGCCCAAAAGGCAAGTTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGY",
        "GCCCGGGGAAGAGGAGGGGCCCGGCGG-AGGCGATAAAAGTGGGGACACATTTCGACGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
        "GGATGGAGAATTAGAGCACTTGTTCTTTTTGCAGAAGCTCAGAATAAACGTTTCGATGGTAA---GGTATTGGCTTACCGTGGCAATGACAGGT-GACGG",
]

###############################################################################
## TESTS: Strict, Sequential

class StrictSequentialDnaTest(extendedtest.ExtendedTestCase):

    def verify(self, dataset):
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqual(len(dataset.tree_lists), 0)
        taxon_set = dataset.taxon_sets[0]
        char_matrix = dataset.char_matrices[0]
        self.assertIs(taxon_set, char_matrix.taxon_set)
        self.assertEqual(len(char_matrix), len(strict_sequential_dna_expected_seq_symbols))
        for i, taxon in enumerate(taxon_set):
            self.assertIn(taxon, char_matrix)
            self.assertIs(char_matrix[taxon], char_matrix[i])
            self.assertEqual(taxon.label, strict_sequential_dna_expected_taxon_labels[i])
            seqs = char_matrix[taxon].symbols_as_string()
            self.assertEqual(seqs, strict_sequential_dna_expected_seq_symbols[i])

    def test_strict_sequential_dna_valid_chars_clean_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        dataset = pr.read(StringIO(strict_sequential_dna_valid_chars_clean_str))
        self.verify(dataset)

    def test_strict_sequential_dna_valid_chars_messy_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        dataset = pr.read(StringIO(strict_sequential_dna_valid_chars_messy_str))
        self.verify(dataset)

    def test_strict_sequential_dna_invalid_chars_messy_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=True)
        dataset = pr.read(StringIO(strict_sequential_dna_invalid_chars_messy_str))
        self.verify(dataset)

    def test_raises_error_strict_sequential_dna_invalid_chars_messy_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=True,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        self.assertRaises(phylip.PhylipReader.PhylipStrictSequentialError,
            pr.read,
            StringIO(strict_sequential_dna_invalid_chars_messy_str))

###############################################################################
## TESTS: Relaxed, Sequential

class RelaxedSequentialDnaTest(extendedtest.ExtendedTestCase):

    def verify(self, dataset, underscores_to_spaces=False):
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqual(len(dataset.tree_lists), 0)
        taxon_set = dataset.taxon_sets[0]
        char_matrix = dataset.char_matrices[0]
        self.assertIs(taxon_set, char_matrix.taxon_set)
        self.assertEqual(len(char_matrix), len(strict_sequential_dna_expected_seq_symbols))
        for i, taxon in enumerate(taxon_set):
            self.assertIn(taxon, char_matrix)
            self.assertIs(char_matrix[taxon], char_matrix[i])
            if underscores_to_spaces:
                expected_label = relaxed_sequential_dna_expected_taxon_labels[i].replace('_', ' ')
            else:
                expected_label = relaxed_sequential_dna_expected_taxon_labels[i]
            self.assertEqual(taxon.label, expected_label)
            seqs = char_matrix[taxon].symbols_as_string()
            self.assertEqual(seqs, relaxed_sequential_dna_expected_seq_symbols[i])

    ### clean, valid ###

    def test_relaxed_single_space_dna_valid_chars_clean_str_with_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_clean_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=False)

    def test_relaxed_single_space_dna_valid_chars_clean_str_no_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_clean_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    def test_relaxed_multi_space_dna_valid_chars_clean_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=True,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_clean_str.replace('$', '  ').replace('_', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    ### messy, valid ###

    def test_relaxed_single_space_dna_valid_chars_messy_str_with_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_messy_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=False)

    def test_relaxed_single_space_dna_valid_chars_messy_str_no_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_messy_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    def test_relaxed_multi_space_dna_valid_chars_messy_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=True,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_sequential_dna_valid_chars_messy_str.replace('$', '  ').replace('_', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    ### messy, invalid ###

    def test_relaxed_single_space_dna_invalid_chars_messy_str_with_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=True)
        s = relaxed_sequential_dna_invalid_chars_messy_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=False)

    def test_relaxed_single_space_dna_invalid_chars_messy_str_no_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=False,
            underscores_to_spaces=True,
            ignore_invalid_chars=True)
        s = relaxed_sequential_dna_invalid_chars_messy_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    def test_relaxed_multi_space_dna_invalid_chars_messy_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=False,
            multispace_delimiter=True,
            underscores_to_spaces=True,
            ignore_invalid_chars=True)
        s = relaxed_sequential_dna_invalid_chars_messy_str.replace('$', '  ').replace('_', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

###############################################################################
## TESTS: Strict, Sequential

class StrictInterleavedDnaTest(extendedtest.ExtendedTestCase):

    def verify(self, dataset):
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqual(len(dataset.tree_lists), 0)
        taxon_set = dataset.taxon_sets[0]
        char_matrix = dataset.char_matrices[0]
        self.assertIs(taxon_set, char_matrix.taxon_set)
        self.assertEqual(len(char_matrix), len(strict_interleaved_dna_expected_seq_symbols))
        for i, taxon in enumerate(taxon_set):
            self.assertIn(taxon, char_matrix)
            self.assertIs(char_matrix[taxon], char_matrix[i])
            self.assertEqual(taxon.label, strict_interleaved_dna_expected_taxon_labels[i])
            seqs = char_matrix[taxon].symbols_as_string()
            self.assertEqual(seqs, strict_interleaved_dna_expected_seq_symbols[i])

    ### clean, valid ###

    def test_strict_interleaved_dna_valid_chars_clean_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=True,
            interleaved=True,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        dataset = pr.read(StringIO(strict_interleaved_dna_valid_chars_clean_str))
        self.verify(dataset)

###############################################################################
## TESTS: Relaxed, Interleaved

class RelaxedInterleavedDnaTest(extendedtest.ExtendedTestCase):

    def verify(self, dataset, underscores_to_spaces=False):
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.char_matrices), 1)
        self.assertEqual(len(dataset.tree_lists), 0)
        taxon_set = dataset.taxon_sets[0]
        char_matrix = dataset.char_matrices[0]
        self.assertIs(taxon_set, char_matrix.taxon_set)
        self.assertEqual(len(char_matrix), len(strict_interleaved_dna_expected_seq_symbols))
        for i, taxon in enumerate(taxon_set):
            self.assertIn(taxon, char_matrix)
            self.assertIs(char_matrix[taxon], char_matrix[i])
            if underscores_to_spaces:
                expected_label = relaxed_interleaved_dna_expected_taxon_labels[i].replace('_', ' ')
            else:
                expected_label = relaxed_interleaved_dna_expected_taxon_labels[i]
            self.assertEqual(taxon.label, expected_label)
            seqs = char_matrix[taxon].symbols_as_string()
            self.assertEqual(seqs, relaxed_interleaved_dna_expected_seq_symbols[i])

    ### clean, valid ###

    def test_relaxed_single_space_dna_valid_chars_clean_str_with_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=True,
            multispace_delimiter=False,
            underscores_to_spaces=False,
            ignore_invalid_chars=False)
        s = relaxed_interleaved_dna_valid_chars_clean_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=False)

    def test_relaxed_single_space_dna_valid_chars_clean_str_no_underscores(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=True,
            multispace_delimiter=False,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_interleaved_dna_valid_chars_clean_str.replace('$', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

    def test_relaxed_multi_space_dna_valid_chars_clean_str(self):
        pr = phylip.PhylipReader(
            data_type='dna',
            strict=False,
            interleaved=True,
            multispace_delimiter=True,
            underscores_to_spaces=True,
            ignore_invalid_chars=False)
        s = relaxed_interleaved_dna_valid_chars_clean_str.replace('$', '  ').replace('_', ' ')
        dataset = pr.read(StringIO(s))
        self.verify(dataset, underscores_to_spaces=True)

class StrictDnaReadWrite(datatest.DataObjectVerificationTestCase):

    def test_round_trip(self):
        ds1 = dendropy.DataSet.get_from_string(strict_sequential_dna_valid_chars_clean_str, 'phylip', data_type='dna', strict=True)
        rw_kwargs = {'strict': True, 'data_type': 'dna'}
        self.roundTripDataSetTest(ds1, 'phylip', writer_kwargs=rw_kwargs, reader_kwargs=rw_kwargs)

class RelaxedReadWrite(datatest.DataObjectVerificationTestCase):

    def test_dna_round_trip(self):
        ds1 = datagen.reference_dna_matrix()
        rw_kwargs = {'data_type': 'dna', 'underscores_to_spaces': True, 'spaces_to_underscores': True}
        self.roundTripData(ds1, 'phylip', writer_kwargs=rw_kwargs, reader_kwargs=rw_kwargs)

    def test_std_round_trip(self):
        ds1 = datagen.reference_standard_binary_matrix()
        rw_kwargs = {'data_type': 'standard', 'underscores_to_spaces': True, 'spaces_to_underscores': True}
        self.roundTripData(ds1, 'phylip', writer_kwargs=rw_kwargs, reader_kwargs=rw_kwargs)

if __name__ == "__main__":
    unittest.main()

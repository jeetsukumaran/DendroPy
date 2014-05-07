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
Tests character sequence map.
"""

import unittest
import dendropy
from dendropy.datamodel import charmatrixmodel
from dendropy.test.support import dendropytest

def get_taxon_namespace(ntax):
    taxon_namespace = dendropy.TaxonNamespace()
    for i in range(ntax):
        label = "T{}".format(i)
        t = taxon_namespace.require_taxon(label=label)
    return taxon_namespace

class CharacterMatrixBasicCRUDTests(dendropytest.ExtendedTestCase):

    def test_setitem_by_taxon(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_idx(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[idx] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_label(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon.label] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_not_in_namespace(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        t = tns[0]
        seq = ["a", "b"]
        with self.assertRaises(ValueError):
            char_matrix[t] = seq
        char_matrix.taxon_namespace.add_taxon(t)
        char_matrix[t] = seq
        self.assertEqual(len(char_matrix), 1)
        self.assertIn(t, char_matrix)
        self.assertEqual(len(char_matrix[t]), len(seq))
        self.assertTrue(isinstance(char_matrix[t], charmatrixmodel.CharacterSequence))
        for c1, c2 in zip(char_matrix[t], seq):
            self.assertEqual(c1, c2)

    def test_setitem_by_idx_not_in_namespace(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(IndexError):
            char_matrix[len(tns)] = []

    def test_setitem_by_idx_not_in_namespace(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(KeyError):
            char_matrix[tns[0].label] = []

    def test_multi_setitem(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        t = tns[0]
        for seq in seqs:
            char_matrix[t] = seq
        for taxon in tns:
            if taxon is t:
                self.assertIn(taxon, char_matrix)
            else:
                self.assertNotIn(taxon, char_matrix)
        seq = seqs[-1]
        self.assertEqual(len(char_matrix), 1)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix[0]), len(seq))
        self.assertTrue(isinstance(char_matrix[0], charmatrixmodel.CharacterSequence))
        for c1, c2 in zip(char_matrix[0], seq):
            self.assertEqual(c1, c2)
        for c1, c2 in zip(char_matrix[0], seqs[1]):
            self.assertNotEqual(c1, c2)

    def test_delitem(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            del char_matrix[taxon]
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
        self.assertEqual(len(char_matrix._taxon_sequence_map), 0)
        self.assertEqual(len(char_matrix), 0)

    def test_clear(self):
        tns = get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        char_matrix.clear()
        self.assertEqual(len(char_matrix._taxon_sequence_map), 0)
        self.assertEqual(len(char_matrix), 0)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)

class CharacterMatrixMetricsTest(dendropytest.ExtendedTestCase):

    def test_sequence_sizes(self):
        seq_sizes = [2, 10, 20, 0, 1]
        tns = get_taxon_namespace(len(seq_sizes))
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), 0)
        self.assertEqual(char_matrix.sequence_size, 0)
        self.assertEqual(char_matrix.max_sequence_size, 0)
        for taxon, seq_size in zip(tns, seq_sizes):
            char_matrix[taxon] = ["x"] * seq_size
        self.assertEqual(len(char_matrix), len(seq_sizes))
        self.assertEqual(char_matrix.sequence_size, seq_sizes[0])
        self.assertEqual(char_matrix.max_sequence_size, max(seq_sizes))

class CharacterMatrixBulkOpsTestCase(dendropytest.ExtendedTestCase):

    def test_fill(self):
        seq_sizes = [2, 10, 20, 0, 1]
        tns = get_taxon_namespace(len(seq_sizes))
        original_sequences = []
        for seq_size in seq_sizes:
            original_sequences.append( ["1"] * seq_size )
        for size in (None, 50, 1, 0, 8):
            for append in (False, True, None):
                kwargs = {}
                if size is None:
                    expected_sizes = [max(seq_sizes)] * len(seq_sizes)
                else:
                    kwargs["size"] = size
                    expected_sizes = [max(size, s) for s in seq_sizes]
                assert len(expected_sizes) == len(original_sequences)
                if append is None:
                    append = True
                else:
                    kwargs["append"] = append
                expected_sequences = []
                for idx, seq in enumerate(original_sequences):
                    if expected_sizes[idx] <= len(seq):
                        expected_sequences.append(list(seq))
                    else:
                        s1 = list(seq)
                        diff = expected_sizes[idx] - len(s1)
                        s2 = ["0"] * diff
                        if append:
                            s = s1 + s2
                        else:
                            s = s2 + s1
                        expected_sequences.append(s)
                    assert len(expected_sequences[idx]) == expected_sizes[idx], \
                            "{}: {}/{}: {}: {} ({})".format(idx, size, append, expected_sequences[idx], len(expected_sequences[idx]), expected_sizes[idx])
                char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
                for taxon, seq in zip(tns, original_sequences):
                    char_matrix[taxon] = seq
                assert len(char_matrix) == len(seq_sizes)
                char_matrix.fill("0", **kwargs)
                for taxon, expected_size, expected_seq in zip(char_matrix, expected_sizes, expected_sequences):
                    obs_seq = char_matrix[taxon]
                    self.assertEqual(len(obs_seq), expected_size)
                    for c1, c2 in zip(obs_seq, expected_seq):
                        self.assertEqual(c1, c2)

if __name__ == "__main__":
    unittest.main()

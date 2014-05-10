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

import collections
import random
import unittest
import dendropy
from dendropy.utility import error
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

class CharacterMatrixFillAndPackTestCase(dendropytest.ExtendedTestCase):

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

    def test_fill_taxa(self):
        tns = get_taxon_namespace(5)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        for taxon in tns[:3]:
            char_matrix[taxon] = "z"
        for taxon in tns[:3]:
            self.assertIn(taxon, char_matrix)
        for taxon in tns[3:]:
            self.assertNotIn(taxon, char_matrix)
        char_matrix.fill_taxa()
        for taxon in tns:
            self.assertIn(taxon, char_matrix)

    def test_fill_taxa(self):
        tns = get_taxon_namespace(5)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        s = ["z"] * 10
        char_matrix[tns[0]] = s
        for taxon in tns[1:3]:
            char_matrix[taxon] = ["x"]
        char_matrix.pack()
        self.assertEqual(len(char_matrix), len(tns))
        for taxon in tns:
            self.assertIn(taxon, char_matrix)
            self.assertEqual(len(char_matrix[taxon]), 10)

class CharacterMatrixBinaryOps(dendropytest.ExtendedTestCase):

    def get_char_matrices(self):
        tns = get_taxon_namespace(3)
        c1 = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        c2 = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        c1[tns[0]] = [1, 1, 1]
        c1[tns[1]] = [2, 2, 2]
        c2[tns[1]] = [3, 3, 3]
        c2[tns[2]] = [4, 4, 4]

        assert len(c1) == 2
        assert tns[0] in c1
        assert tns[1] in c1
        assert tns[2] not in c1

        assert len(c2) == 2
        assert tns[0] not in c2
        assert tns[1] in c2
        assert tns[2] in c2

        return c1, c2, tns

    def verify_sequence_equal(self, s1, s2, expected_length=None):
        if expected_length is not None:
            self.assertEqual(len(s1), expected_length)
            self.assertEqual(len(s2), expected_length)
        self.assertEqual(len(s1), len(s2))
        self.assertIsNot(s1, s2)
        for c1, c2 in zip(s1, s2):
            self.assertEqual(c1, c2)

    def verify_independent_matrices(self, c1, c2):
        assert c1.taxon_namespace is c2.taxon_namespace
        for taxon in c1.taxon_namespace:
            if taxon in c1 and taxon in c2:
                self.assertIsNot(c1[taxon], c2[taxon])

    def test_add_sequences_fail(self):
        c1 = charmatrixmodel.CharacterMatrix()
        c2 = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(error.TaxonNamespaceError):
            c1.add_sequences(c2)

    def test_add_sequences(self):
        c1, c2, tns = self.get_char_matrices()
        c1.add_sequences(c2)
        self.verify_independent_matrices(c1, c2)
        self.assertEqual(len(c1), 3)
        self.assertIn(tns[0], c1)
        self.assertIn(tns[1], c1)
        self.assertIn(tns[2], c1)
        self.verify_sequence_equal(c1[tns[0]], [1, 1, 1])
        self.verify_sequence_equal(c1[tns[1]], [2, 2, 2])
        self.verify_sequence_equal(c1[tns[2]], [4, 4, 4])

    def test_replace_sequences_fail(self):
        c1 = charmatrixmodel.CharacterMatrix()
        c2 = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(error.TaxonNamespaceError):
            c1.replace_sequences(c2)

    def test_replace_sequences(self):
        c1, c2, tns = self.get_char_matrices()
        c1.replace_sequences(c2)
        self.verify_independent_matrices(c1, c2)
        self.assertEqual(len(c1), 2)
        self.assertIn(tns[0], c1)
        self.assertIn(tns[1], c1)
        self.assertNotIn(tns[2], c1)
        self.verify_sequence_equal(c1[tns[0]], [1, 1, 1])
        self.verify_sequence_equal(c1[tns[1]], [3, 3, 3])

    def test_update_sequences_fail(self):
        c1 = charmatrixmodel.CharacterMatrix()
        c2 = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(error.TaxonNamespaceError):
            c1.update_sequences(c2)

    def test_update_sequences(self):
        c1, c2, tns = self.get_char_matrices()
        c1.update_sequences(c2)
        self.verify_independent_matrices(c1, c2)
        self.assertEqual(len(c1), 3)
        self.assertIn(tns[0], c1)
        self.assertIn(tns[1], c1)
        self.assertIn(tns[2], c1)
        self.verify_sequence_equal(c1[tns[0]], [1, 1, 1])
        self.verify_sequence_equal(c1[tns[1]], [3, 3, 3])
        self.verify_sequence_equal(c1[tns[2]], [4, 4, 4])

    def test_extend_sequences_fail(self):
        c1 = charmatrixmodel.CharacterMatrix()
        c2 = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(error.TaxonNamespaceError):
            c1.extend_sequences(c2)

    def test_extend_sequences(self):
        c1, c2, tns = self.get_char_matrices()
        c1.extend_sequences(c2)
        self.verify_independent_matrices(c1, c2)
        self.assertEqual(len(c1), 2)
        self.assertIn(tns[0], c1)
        self.assertIn(tns[1], c1)
        self.assertNotIn(tns[2], c1)
        self.verify_sequence_equal(c1[tns[0]], [1, 1, 1])
        self.verify_sequence_equal(c1[tns[1]], [2, 2, 2, 3, 3, 3])

    def test_extend_matrix_fail(self):
        c1 = charmatrixmodel.CharacterMatrix()
        c2 = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(error.TaxonNamespaceError):
            c1.extend_matrix(c2)

    def test_extend_matrix(self):
        c1, c2, tns = self.get_char_matrices()
        c1.extend_matrix(c2)
        self.verify_independent_matrices(c1, c2)
        self.assertEqual(len(c1), 3)
        self.assertIn(tns[0], c1)
        self.assertIn(tns[1], c1)
        self.assertIn(tns[2], c1)
        self.verify_sequence_equal(c1[tns[0]], [1, 1, 1])
        self.verify_sequence_equal(c1[tns[1]], [2, 2, 2, 3, 3, 3])
        self.verify_sequence_equal(c1[tns[2]], [4, 4, 4])

class CharacterMatrixTaxonManagement(dendropytest.ExtendedTestCase):

    def test_assign_taxon_namespace(self):
        tns = get_taxon_namespace(5)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertIs(char_matrix.taxon_namespace, tns)

class CharacterMatrixIteratorTests(dendropytest.ExtendedTestCase):

    def setUp(self):
        self.rng = random.Random()

    def test_standard_iterator(self):
        tns = get_taxon_namespace(100)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        taxa = list(tns)
        self.rng.shuffle(taxa)
        included = set()
        excluded = set()
        for idx, taxon in enumerate(taxa):
            if self.rng.uniform(0, 1) < 0.5:
                included.add(taxon)
                char_matrix[taxon] = [0]
            else:
                excluded.add(taxon)
        expected = [taxon for taxon in tns if taxon in included]
        self.assertEqual(len(char_matrix), len(expected))
        observed = [taxon for taxon in char_matrix]
        self.assertEqual(observed, expected)

class CharacterMatrixIdentity(unittest.TestCase):

    def setUp(self):
        self.tns = dendropy.TaxonNamespace()
        self.t1 = charmatrixmodel.CharacterMatrix(label="a", taxon_namespace=self.tns)
        self.t2 = charmatrixmodel.CharacterMatrix(label="a", taxon_namespace=self.tns)
        self.t3 = charmatrixmodel.CharacterMatrix(label="a")

    def test_equal(self):
        self.assertNotEqual(self.t1, self.t2)

    def test_hash_dict_membership(self):
        k = {}
        k[self.t1] = 1
        k[self.t2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.t1], 1)
        self.assertEqual(k[self.t2], 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        del k[self.t1]
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.t1)
        k.add(self.t2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        k.discard(self.t1)
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

class TestCharacterMatrixUpdateTaxonNamespace(
        dendropytest.ExtendedTestCase):

    def get_char_matrix(self):
        labels = [
            "z01" , "<NONE>" , "z03" , "z04" , "z05" ,
            "z06" , None     , None  , "z09" , "z10" ,
            "z11" , "<NONE>" , None  , "z14" , "z15" ,
                ]
        char_matrix = charmatrixmodel.CharacterMatrix()
        char_matrix.expected_labels = []
        char_matrix.expected_taxa = set()
        for label in labels:
            t = dendropy.Taxon(label=None)
            char_matrix.taxon_namespace.add_taxon(t)
            char_matrix[t] = [1,1,1]
            char_matrix.expected_taxa.add(t)
            char_matrix.expected_labels.append(t.label)
        char_matrix.taxon_namespace = dendropy.TaxonNamespace()
        return char_matrix

    def test_update(self):
        char_matrix = self.get_char_matrix()
        char_matrix.taxon_namespace = dendropy.TaxonNamespace()
        original_tns = char_matrix.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        char_matrix.update_taxon_namespace()
        char_matrix.update_taxon_namespace()
        char_matrix.update_taxon_namespace()
        self.assertIs(char_matrix.taxon_namespace, original_tns)
        self.assertEqual(len(char_matrix.taxon_namespace), len(char_matrix.expected_taxa))
        for taxon in char_matrix:
            self.assertIn(taxon, char_matrix.taxon_namespace)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        self.assertCountEqual(new_taxa, char_matrix.expected_taxa)
        self.assertCountEqual(new_labels, char_matrix.expected_labels)

class TestCharacterMatrixReconstructAndMigrateTaxonNamespace(
        dendropytest.ExtendedTestCase):

    def get_char_matrix(self, labels=None):
        char_matrix = charmatrixmodel.CharacterMatrix()
        if labels is None:
            labels = [str(i) for i in range(1000)]
        char_matrix.expected_labels = []
        char_matrix.original_taxa = []
        self.rng.shuffle(labels)
        for label in labels:
            t = dendropy.Taxon(label=label)
            char_matrix.taxon_namespace.add_taxon(t)
            char_matrix.original_taxa.append(t)
            char_matrix[t].original_taxon = t
            char_matrix.expected_labels.append(label)
            seq = [self.rng.randint(0, 100) for _ in range(4)]
            char_matrix[t] = seq
            char_matrix[t].original_seq = char_matrix[t]
            char_matrix[t].original_taxon = t
            char_matrix[t].label = label
        assert len(char_matrix.taxon_namespace) == len(char_matrix.original_taxa)
        return char_matrix

    def setUp(self):
        self.rng = random.Random()
        labels = [
                "a", "a", "2", "2", "b", "B",
                "B", "h", "H", "h", None, None,
                "H", "J", "j",
                ]
        self.char_matrix = self.get_char_matrix(labels=labels)

    def verify_taxon_namespace_reconstruction(self,
            char_matrix=None,
            unify_taxa_by_label=False,
            case_sensitive_label_mapping=True,
            original_tns=None):
        if char_matrix is None:
            char_matrix = self.char_matrix
        if unify_taxa_by_label:
            if not case_sensitive_label_mapping:
                expected_labels = list(set((label.upper() if label is not None else None) for label in char_matrix.expected_labels))
            else:
                expected_labels = list(set(label for label in char_matrix.expected_labels))
        else:
            expected_labels = [label for label in char_matrix.expected_labels]
        seen_taxa = []
        for taxon in char_matrix:
            seq = char_matrix[taxon]
            self.assertIsNot(taxon, seq.original_taxon)
            if not case_sensitive_label_mapping and taxon.label is not None:
                self.assertEqual(taxon.label.upper(), seq.original_taxon.label.upper())
                self.assertEqual(seq.label.upper(), taxon.label.upper())
            else:
                self.assertEqual(taxon.label, seq.original_taxon.label)
                self.assertEqual(seq.label, taxon.label)
            self.assertNotIn(seq.original_taxon, char_matrix.taxon_namespace)
            self.assertIn(seq.original_taxon, char_matrix.original_taxa)
            self.assertIn(taxon, char_matrix.taxon_namespace)
            self.assertNotIn(taxon, char_matrix.original_taxa)
            if original_tns is not None:
                self.assertNotIn(taxon, original_tns)
            if taxon not in seen_taxa:
                seen_taxa.append(taxon)
            else:
                self.assertTrue(unify_taxa_by_label)
                if not case_sensitive_label_mapping:
                    self.assertIn(taxon.label, [t.label for t in seen_taxa])
                else:
                    if taxon.label is None:
                        self.assertIs(seq.original_taxon.label, None)
                        self.assertEqual([t.label for t in seen_taxa].count(None), 1)
                    else:
                        x1 = [t.label.upper() for t in seen_taxa if t.label is not None]
                        self.assertIn(taxon.label.upper(), x1)
        self.assertEqual(len(seen_taxa), len(char_matrix.taxon_namespace))
        if not case_sensitive_label_mapping:
            seen_labels = [(t.label.upper() if t.label is not None else None) for t in seen_taxa]
        else:
            seen_labels = [t.label for t in seen_taxa]
        c1 = collections.Counter(expected_labels)
        c2 = collections.Counter(seen_labels)
        self.assertEqual(c2-c1, collections.Counter())
        self.assertEqual(c1-c2, collections.Counter())
        self.assertEqual(c1, c2)
        self.assertEqual(len(char_matrix.taxon_namespace), len(expected_labels))
        if not unify_taxa_by_label:
            self.assertEqual(len(char_matrix.taxon_namespace), len(char_matrix.original_taxa))

    def test_basic_reconstruction(self):
        char_matrix = self.get_char_matrix()
        tns = char_matrix.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        char_matrix.taxon_namespace = new_tns
        char_matrix.reconstruct_taxon_namespace(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True)
        # print("\n--\n")
        # for t in self.char_matrix:
            # print("{}: {}".format(repr(t), self.char_matrix[t]))
            # assert t in self.char_matrix.taxon_namespace
        self.assertIsNot(char_matrix.taxon_namespace, tns)
        self.assertIs(char_matrix.taxon_namespace, new_tns)
        if len(char_matrix.taxon_namespace) != len(tns):
            x1 = [t.label for t in char_matrix.taxon_namespace]
            x2 = [t.label for t in tns]
            c1 = collections.Counter(x1)
            c2 = collections.Counter(x2)
            c3 = c2 - c1
            print(c3)
        self.assertEqual(len(char_matrix.taxon_namespace), len(tns))
        original_labels = [t.label for t in tns]
        new_labels = [t.label for t in new_tns]
        self.assertCountEqual(new_labels, original_labels)
        for taxon in char_matrix:
            self.assertIs(char_matrix[taxon], char_matrix[taxon].original_seq)
            self.assertIn(taxon, char_matrix.taxon_namespace)
            self.assertNotIn(taxon, tns)

    def test_reconstruct_taxon_namespace_non_unifying(self):
        original_tns = self.char_matrix.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        self.char_matrix._taxon_namespace = new_tns
        self.assertEqual(len(self.char_matrix.taxon_namespace), 0)
        self.char_matrix.reconstruct_taxon_namespace(unify_taxa_by_label=False,
                case_sensitive_label_mapping=True)
        self.assertIsNot(self.char_matrix.taxon_namespace, original_tns)
        self.assertIs(self.char_matrix.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True)

#     def test_reconstruct_taxon_namespace_unifying_case_sensitive(self):
#         original_tns = self.tree_list.taxon_namespace
#         new_tns = dendropy.TaxonNamespace()
#         self.tree_list._taxon_namespace = new_tns
#         self.assertEqual(len(self.tree_list.taxon_namespace), 0)
#         self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=True)
#         self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
#         self.assertIs(self.tree_list.taxon_namespace, new_tns)
#         self.verify_taxon_namespace_reconstruction(
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=True,
#                 original_tns=original_tns)

#     def test_reconstruct_taxon_namespace_unifying_case_insensitive(self):
#         original_tns = self.tree_list.taxon_namespace
#         new_tns = dendropy.TaxonNamespace()
#         self.tree_list._taxon_namespace = new_tns
#         self.assertEqual(len(self.tree_list.taxon_namespace), 0)
#         self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=False)
#         self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
#         self.assertIs(self.tree_list.taxon_namespace, new_tns)
#         self.verify_taxon_namespace_reconstruction(
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=False,
#                 original_tns=original_tns)

#     def test_basic_migration(self):
#         tns = dendropy.TaxonNamespace()
#         trees = []
#         for idx in range(5):
#             tree, anodes, lnodes, inodes = self.get_tree(
#                     suppress_internal_node_taxa=False,
#                     suppress_leaf_node_taxa=False,
#                     taxon_namespace=tns)
#             trees.append(tree)
#         tree_list = dendropy.TreeList(taxon_namespace=tns)
#         tree_list._trees = trees
#         new_tns = dendropy.TaxonNamespace()
#         tree_list.taxon_namespace = new_tns
#         tree_list.migrate_taxon_namespace(
#                 new_tns,
#                 unify_taxa_by_label=False,
#                 case_sensitive_label_mapping=True)
#         self.assertIsNot(tree_list.taxon_namespace, tns)
#         self.assertIs(tree_list.taxon_namespace, new_tns)
#         self.assertEqual(len(tree_list.taxon_namespace), len(tns))
#         original_labels = [t.label for t in tns]
#         new_labels = [t.label for t in new_tns]
#         self.assertCountEqual(new_labels, original_labels)
#         for tree in tree_list:
#             self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
#             for nd in tree:
#                 if nd.taxon is not None:
#                     self.assertIn(nd.taxon, tree.taxon_namespace)
#                     self.assertNotIn(nd.taxon, tns)

#     def test_migrate_taxon_namespace_non_unifying(self):
#         original_tns = self.tree_list.taxon_namespace
#         new_tns = dendropy.TaxonNamespace()
#         self.tree_list.migrate_taxon_namespace(
#                 new_tns,
#                 unify_taxa_by_label=False,
#                 case_sensitive_label_mapping=True)
#         self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
#         self.assertIs(self.tree_list.taxon_namespace, new_tns)
#         self.verify_taxon_namespace_reconstruction(
#                 unify_taxa_by_label=False,
#                 case_sensitive_label_mapping=True,
#                 original_tns=original_tns)

#     def test_migrate_taxon_namespace_unifying_case_sensitive(self):
#         original_tns = self.tree_list.taxon_namespace
#         new_tns = dendropy.TaxonNamespace()
#         self.tree_list.migrate_taxon_namespace(
#                 new_tns,
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=True)
#         self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
#         self.assertIs(self.tree_list.taxon_namespace, new_tns)
#         self.verify_taxon_namespace_reconstruction(
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=True,
#                 original_tns=original_tns)

#     def test_migrate_taxon_namespace_unifying_case_insensitive(self):
#         original_tns = self.tree_list.taxon_namespace
#         new_tns = dendropy.TaxonNamespace()
#         self.tree_list.migrate_taxon_namespace(
#                 new_tns,
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=False)
#         self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
#         self.assertIs(self.tree_list.taxon_namespace, new_tns)
#         self.verify_taxon_namespace_reconstruction(
#                 unify_taxa_by_label=True,
#                 case_sensitive_label_mapping=False,
#                 original_tns=original_tns)

# class TestTreeListAppend(
#         datagen_curated_test_tree.CuratedTestTree,
#         unittest.TestCase):

#     def setUp(self):
#         self.native_tns = dendropy.TaxonNamespace()
#         self.tree_list = dendropy.TreeList(taxon_namespace=self.native_tns)
#         self.foreign_tns = dendropy.TaxonNamespace()
#         self.foreign_tree, anodes, lnodes, inodes = self.get_tree(
#                 suppress_internal_node_taxa=False,
#                 suppress_leaf_node_taxa=False,
#                 taxon_namespace=self.foreign_tns)
#         for nd in self.foreign_tree:
#             nd.original_taxon = nd.taxon
#         self.check_tns = dendropy.TaxonNamespace()
#         self.check_tree, anodes, lnodes, inodes = self.get_tree(
#                 suppress_internal_node_taxa=False,
#                 suppress_leaf_node_taxa=False,
#                 taxon_namespace=self.check_tns)

#     def test_append_default(self):
#         self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
#         self.tree_list.append(self.foreign_tree)
#         self.assertEqual(len(self.tree_list), 1)
#         self.assertIn(self.foreign_tree, self.tree_list)
#         self.assertIs(self.foreign_tree, self.tree_list[0])
#         self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
#         self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
#         self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
#         for nd in self.foreign_tree:
#             if nd.taxon:
#                 self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
#                 self.assertIsNot(nd.taxon, nd.original_taxon)
#                 self.assertIn(nd.original_taxon, self.foreign_tns)
#                 self.assertNotIn(nd.original_taxon, self.tree_list.taxon_namespace)
#                 self.assertEqual(nd.taxon.label, nd.original_taxon.label)

#     def test_append_migrate_matching_labels(self):
#         kwargs_groups = [
#                 {"taxon_import_strategy": "migrate", "unify_taxa_by_label": True},
#                 {"taxon_import_strategy": "migrate", "unify_taxa_by_label": False},
#                 {"taxon_import_strategy": "add", },
#         ]
#         for kwargs in kwargs_groups:
#             self.setUp()
#             self.assertEqual(len(self.tree_list.taxon_namespace), 0)
#             native_tree, anodes, lnodes, inodes = self.get_tree(
#                     suppress_internal_node_taxa=False,
#                     suppress_leaf_node_taxa=False,
#                     taxon_namespace=self.native_tns)
#             self.assertEqual(len(self.tree_list.taxon_namespace), len(self.postorder_sequence))
#             self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
#             original_tns_len = len(self.tree_list.taxon_namespace)
#             self.tree_list.append(self.foreign_tree, **kwargs)
#             self.assertEqual(len(self.tree_list), 1)
#             self.assertIn(self.foreign_tree, self.tree_list)
#             self.assertIs(self.foreign_tree, self.tree_list[0])
#             self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
#             if kwargs["taxon_import_strategy"] == "add":
#                 self.assertEqual(len(self.tree_list.taxon_namespace),
#                         original_tns_len + len(self.foreign_tns))
#                 for nd in self.foreign_tree:
#                     self.assertIn(nd.taxon, self.foreign_tns)
#                     self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
#             else:
#                 if "unify_taxa_by_label" not in kwargs or not kwargs["unify_taxa_by_label"]:
#                     self.assertEqual(len(self.tree_list.taxon_namespace),
#                             original_tns_len + len(self.foreign_tns))
#                 else:
#                     self.assertEqual(len(self.tree_list.taxon_namespace), original_tns_len)
#                 for nd in self.foreign_tree:
#                     self.assertNotIn(nd.taxon, self.foreign_tns)
#                     self.assertIn(nd.taxon, self.tree_list.taxon_namespace)

#     def test_append_add(self):
#         self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
#         self.tree_list.append(self.foreign_tree,
#                 taxon_import_strategy="add")
#         self.assertEqual(len(self.tree_list), 1)
#         self.assertIn(self.foreign_tree, self.tree_list)
#         self.assertIs(self.foreign_tree, self.tree_list[0])
#         self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
#         self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
#         self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
#         for nd in self.foreign_tree:
#             if nd.taxon:
#                 self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
#                 self.assertIs(nd.taxon, nd.original_taxon)
#                 self.assertIn(nd.original_taxon, self.foreign_tns)
#                 self.assertIn(nd.original_taxon, self.tree_list.taxon_namespace)

# class TreeListCreation(unittest.TestCase):

#     def test_create_with_taxon_namespace(self):
#         tns = dendropy.TaxonNamespace()
#         tt = TreeList(label="a", taxon_namespace=tns)
#         self.assertEqual(tt.label, "a")
#         self.assertIs(tt.taxon_namespace, tns)

# class TreeListCreatingAndCloning(
#         datagen_curated_test_tree.CuratedTestTree,
#         compare_and_validate.Comparator,
#         unittest.TestCase):

#     def add_tree_annotations(self, tree):
#         for idx, nd in enumerate(tree):
#             if idx % 2 == 0:
#                 nd.edge.label = "E{}".format(idx)
#                 nd.edge.length = idx
#             an1 = nd.annotations.add_new("a{}".format(idx),
#                     "{}{}{}".format(nd.label, nd.taxon, idx))
#             an2 = nd.annotations.add_bound_attribute("label")
#             an3 = an1.annotations.add_bound_attribute("name")
#             ae1 = nd.edge.annotations.add_new("a{}".format(idx),
#                     "{}{}".format(nd.edge.label, idx))
#             ae2 = nd.edge.annotations.add_bound_attribute("label")
#             ae3 = ae1.annotations.add_bound_attribute("name")
#         tree.annotations.add_new("a", 0)
#         tree.label = "hello"
#         b = tree.annotations.add_bound_attribute("label")
#         b.annotations.add_new("c", 3)

#     def add_tree_list_annotations(self, tree_list):
#         tree_list.annotations.add_new("a", 0)
#         tree_list.label = "hello"
#         b = tree_list.annotations.add_bound_attribute("label")
#         b.annotations.add_new("c", 3)

#     def add_taxon_namespace_annotations(self, tns):
#         for idx, taxon in enumerate(tns):
#             a = taxon.annotations.add_new("!color", str(idx))
#             a.annotations.add_new("setbytest", "a")

#     def setUp(self):
#         self.num_trees = 5
#         tree1, anodes1, lnodes1, inodes1 = self.get_tree(
#                 suppress_internal_node_taxa=False,
#                 suppress_leaf_node_taxa=False)
#         self.original_taxon_labels = [t.label for t in tree1.taxon_namespace]
#         assert len(self.original_taxon_labels) == len(anodes1)

#     def get_tree_list(self):
#         tlist1 = TreeList()
#         self.num_trees = 5
#         for idx in range(self.num_trees):
#             tree1, anodes1, lnodes1, inodes1 = self.get_tree(
#                     suppress_internal_node_taxa=False,
#                     suppress_leaf_node_taxa=False,
#                     taxon_namespace=tlist1.taxon_namespace)
#             self.add_tree_annotations(tree1)
#             tlist1.append(tree1)
#         self.add_tree_list_annotations(tlist1)
#         self.add_taxon_namespace_annotations(tlist1.taxon_namespace)
#         return tlist1

#     def test_shallow_copy_with_initializer_list(self):
#         tlist1 = self.get_tree_list()
#         trees = tlist1._trees
#         tlist2 = dendropy.TreeList(trees)
#         self.assertEqual(len(tlist2), self.num_trees)
#         for tcopy, toriginal in zip(tlist2, trees):
#             self.assertIs(tcopy, toriginal)
#             self.assertIs(tcopy.taxon_namespace, tlist2.taxon_namespace)

#     def test_clone0(self):
#         tlist1 = self.get_tree_list()
#         for tlist2 in (
#                 tlist1.clone(0),
#                 ):
#             self.assertIs(tlist2.taxon_namespace, tlist1.taxon_namespace)
#             self.assertEqual(len(tlist2), self.num_trees)
#             for tcopy, toriginal in zip(tlist2, tlist1):
#                 self.assertIs(tcopy, toriginal)
#                 self.assertIs(tcopy.taxon_namespace, tlist2.taxon_namespace)

#     def test_taxon_namespace_scoped_copy(self):
#         tlist1 = self.get_tree_list()
#         for tlist2 in (
#                 tlist1.clone(1),
#                 dendropy.TreeList(tlist1),
#                 tlist1.taxon_namespace_scoped_copy(),):
#             self.compare_distinct_tree_list(tlist2, tlist1,
#                     taxon_namespace_scoped=True,
#                     compare_tree_annotations=True,
#                     compare_taxon_annotations=True)

#     def test_deepcopy_including_namespace(self):
#         tlist1 = self.get_tree_list()
#         for idx, tlist2 in enumerate((
#                 tlist1.clone(2),
#                 copy.deepcopy(tlist1),
#                 )):
#             self.compare_distinct_tree_list(tlist2, tlist1,
#                     taxon_namespace_scoped=False,
#                     compare_tree_annotations=True,
#                     compare_taxon_annotations=True)

#     def test_deepcopy_excluding_namespace(self):
#         tlist1 = self.get_tree_list()
#         tlist2 = dendropy.TreeList(tlist1,
#                 taxon_namespace=dendropy.TaxonNamespace())
#         self.compare_distinct_tree_list(tlist2, tlist1,
#                 taxon_namespace_scoped=False,
#                 compare_tree_annotations=True,
#                 compare_taxon_annotations=False)

if __name__ == "__main__":
    unittest.main()

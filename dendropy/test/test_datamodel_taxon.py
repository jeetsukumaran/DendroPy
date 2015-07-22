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
Tests creation, reading, update, deletion of Taxon and TaxonNamespace objects.
"""

import collections
import unittest
import copy
from dendropy import Taxon, TaxonNamespace
from dendropy.test.support import compare_and_validate

class TaxonIdentity(compare_and_validate.Comparator, unittest.TestCase):

    def setUp(self):
        self.t1 = Taxon("a")
        self.t2 = Taxon("a")

    def test_equal(self):
        # two distinct |Taxon| objects are never equal, even if all
        # member values are the same.
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

# note that compare_and_validate.Comparator must be listed first,
# otherwise setUp will not be called
class TaxonCloning(compare_and_validate.Comparator, unittest.TestCase):

    def test_construct_from_another(self):
        t1 = Taxon("a")
        for t2 in (Taxon(t1), copy.deepcopy(t1), t1.clone(2)):
            self.assertIsNot(t1, t2)
            self.assertNotEqual(t1, t2)
            self.assertEqual(t1.label, t2.label)

    def test_construct_from_another_with_simple_annotations(self):
        t1 = Taxon("a")
        t1.annotations.add_new("a", 0)
        t1.annotations.add_new("b", 1)
        t1.annotations.add_new("c", 3)
        for t2 in (Taxon(t1), copy.deepcopy(t1), t1.clone(2)):
            self.assertIsNot(t1, t2)
            self.assertNotEqual(t1, t2)
            self.assertEqual(t1.label, t2.label)
            self.assertTrue(hasattr(t1, "annotations"))
            self.assertTrue(hasattr(t2, "annotations"))
            self.assertEqual(len(t1.annotations), len(t2.annotations))
            self.compare_distinct_annotables(t1, t2)

    def test_construct_from_another_with_complex_annotations(self):
        t1 = Taxon("a")
        t1.annotations.add_new("a", 0)
        b = t1.annotations.add_new("b", (t1, "label"), is_attribute=True)
        b.annotations.add_new("c", 3)
        for t2 in (Taxon(t1), copy.deepcopy(t1), t1.clone(2)):
            self.assertIsNot(t1, t2)
            self.assertNotEqual(t1, t2)
            self.assertEqual(t1.label, t2.label)
            self.assertTrue(hasattr(t1, "annotations"))
            self.assertTrue(hasattr(t2, "annotations"))
            self.assertEqual(len(t1.annotations), len(t2.annotations))
            self.compare_distinct_annotables(t1, t2)
            t1.label = "x"
            self.assertEqual(t1.annotations[1].value, "x")
            self.assertEqual(t2.annotations[1].value, "a")
            t2.label = "z"
            self.assertEqual(t1.annotations[1].value, "x")
            self.assertEqual(t2.annotations[1].value, "z")
            t1.label = "a"

    def test_simple_copy(self):
        t1 = Taxon("a")
        with self.assertRaises(TypeError):
            copy.copy(t1)
        with self.assertRaises(TypeError):
            t1.clone(0)

    def test_taxon_namespace_scoped_copy(self):
        t1 = Taxon("a")
        for t2 in (t1.clone(1), t1.taxon_namespace_scoped_copy()):
            self.assertIs(t2, t1)

class TaxonNamespaceTaxonManagement(unittest.TestCase):

    def setUp(self):
        self.str_labels = ["a", "a", "b", "c", "d", "e", "_", "_", "_", "z", "z", "z"]
        self.taxa = [
                Taxon("t1"), Taxon("t2"), Taxon("t3"),
                ]
        self.taxa_labels = [t.label for t in self.taxa]

    def validate_taxon_concepts(self, tns, labels, respect_order=False):
        expected_labels = labels
        discovered_labels = []
        for t in tns._taxa:
            discovered_labels.append(t.label)
        self.assertEqual(len(discovered_labels), len(expected_labels))
        if respect_order:
            for x1, x2 in zip(discovered_labels, expected_labels):
                self.assertEqual(x1, x2)
        else:
            c1 = collections.Counter(discovered_labels)
            c2 = collections.Counter(expected_labels)
            self.assertEqual(c1, c2)

    ### initialization ###

    def test_initialize_from_str_list(self):
        tns = TaxonNamespace(self.str_labels)
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_initialize_from_taxon_list(self):
        tns = TaxonNamespace(self.taxa)
        self.validate_taxon_concepts(tns, self.taxa_labels)
        for t in self.taxa:
            self.assertIn(t, tns._taxa)

    def test_initialize_other_taxon_namespace(self):
        tns1 = TaxonNamespace(self.taxa)
        tns2 = TaxonNamespace(tns1)
        self.assertIsNot(tns1, tns2)
        self.validate_taxon_concepts(tns1, self.taxa_labels)
        self.validate_taxon_concepts(tns2, self.taxa_labels)
        for t in self.taxa:
            self.assertIn(t, tns1._taxa)
            self.assertIn(t, tns2._taxa)
        for t1, t2 in zip(tns1, tns2):
            self.assertIs(t1, t2)
        # self.assertEqual(tns1, tns2)

    ### adding ###

    def test_basic_adding(self):
        tns = TaxonNamespace()
        self.assertEqual(len(tns), 0)
        for idx, label in enumerate(self.str_labels):
            tns.add_taxon(Taxon(label=label))
            self.assertEqual(len(tns), idx+1)
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_basic_adding_to_immutable(self):
        tns = TaxonNamespace()
        self.assertEqual(len(tns), 0)
        tns.is_mutable = False
        for idx, label in enumerate(self.str_labels):
            with self.assertRaises(TypeError):
                tns.add_taxon(Taxon(label=label))
            self.assertEqual(len(tns), 0)

    def test_add_taxon(self):
        tns = TaxonNamespace()
        for t in self.taxa:
            tns.add_taxon(t)
        self.validate_taxon_concepts(tns, self.taxa_labels)
        for t in self.taxa:
            self.assertIn(t, tns._taxa)

    def test_add_taxa(self):
        tns = TaxonNamespace()
        tns.add_taxa(self.taxa)
        self.validate_taxon_concepts(tns, self.taxa_labels)
        for t in self.taxa:
            self.assertIn(t, tns._taxa)

    def test_add_taxon_duplicate(self):
        tns = TaxonNamespace(self.taxa)
        self.validate_taxon_concepts(tns, self.taxa_labels)
        tns.add_taxon(self.taxa[0])
        self.assertEqual(len(tns), len(self.taxa))
        for t1, t2 in zip(tns, self.taxa):
            self.assertIs(t1, t2)

    def test_add_taxa_duplicate(self):
        tns = TaxonNamespace(self.taxa)
        self.validate_taxon_concepts(tns, self.taxa_labels)
        tns.add_taxa(self.taxa)
        self.assertEqual(len(tns), len(self.taxa))
        for t1, t2 in zip(tns, self.taxa):
            self.assertIs(t1, t2)

    def test_new_taxon(self):
        tns = TaxonNamespace()
        for idx, label in enumerate(self.str_labels):
            t = tns.new_taxon(label)
            self.assertTrue(isinstance(t, Taxon))
            self.assertEqual(t.label, label)
            self.assertEqual(len(tns), idx+1)
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_new_taxa(self):
        tns = TaxonNamespace()
        tns.new_taxa(self.str_labels)
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_new_taxon_to_immutable(self):
        tns = TaxonNamespace()
        tns.is_mutable = False
        for idx, label in enumerate(self.str_labels):
            with self.assertRaises(TypeError):
                t = tns.new_taxon(label)
            self.assertEqual(len(tns), 0)

    def test_new_taxa_to_immutable(self):
        tns = TaxonNamespace()
        tns.is_mutable = False
        with self.assertRaises(TypeError):
            tns.new_taxa(self.str_labels)
        self.assertEqual(len(tns), 0)

    ### access ###

    def test_len(self):
        tns = TaxonNamespace(self.str_labels)
        self.assertEqual(len(tns), len(self.str_labels))

    def test_getitem_by_index(self):
        tns = TaxonNamespace(self.str_labels)
        for index, label in enumerate(self.str_labels):
            t = tns[index]
            self.assertEqual(t.label, label)

    def test_getitem_by_index_error(self):
        tns = TaxonNamespace(self.str_labels)
        with self.assertRaises(IndexError):
            t = tns[len(self.str_labels)+1]

    def test_getitem_slices(self):
        tns = TaxonNamespace(self.str_labels)
        slices = [
                (-1,1,-1),
                (2,4),
                (1,),
                (None,1)
                ]
        for s in slices:
            labels = self.str_labels[slice(*s)]
            taxa = tns[slice(*s)]
            taxa_labels = [t.label for t in taxa]
            self.assertEqual(taxa_labels, labels)

    def test_getitem_by_label_error(self):
        tns = TaxonNamespace(self.str_labels)
        with self.assertRaises(ValueError):
            tns[self.str_labels[0]]

#     def test_getitem_by_label_error(self):
#         tns = TaxonNamespace(self.str_labels)
#         check = ["u", "x", "y",]
#         for label in check:
#             assert label not in self.str_labels
#             with self.assertRaises(LookupError):
#                 t = tns[label]

    def test_contains_taxa(self):
        tns = TaxonNamespace(self.taxa)
        for taxon in self.taxa:
            self.assertIn(taxon, tns)

    def test_no_contains_taxa(self):
        tns = TaxonNamespace(self.taxa)
        taxa2 = [Taxon(label=t.label) for t in self.taxa]
        for taxon in taxa2:
            self.assertNotIn(taxon, tns)

    def test_has_label(self):
        tns = TaxonNamespace(self.str_labels)
        for label in self.str_labels:
            self.assertTrue(tns.has_taxon_label(label))

    def test_no_has_label(self):
        tns = TaxonNamespace(self.str_labels)
        check = ["u", "x", "y",]
        for label in check:
            assert label not in self.str_labels
            self.assertFalse(tns.has_taxon_label(label))

    def test_has_label_case_sensitivity(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        assert labels_upper
        for label in labels_upper:
            tns.is_case_sensitive = True
            self.assertFalse(tns.has_taxon_label(label))
            tns.is_case_sensitive = False
            self.assertTrue(tns.has_taxon_label(label))

    def test_has_labels(self):
        tns = TaxonNamespace(self.str_labels)
        self.assertTrue(tns.has_taxa_labels(self.str_labels))

    def test_no_has_labels(self):
        tns = TaxonNamespace(self.str_labels)
        check = ["u", "x", "y",]
        for label in check:
            assert label not in self.str_labels
        self.assertFalse(tns.has_taxa_labels(check))
        self.assertFalse(tns.has_taxa_labels(check + self.str_labels))

    def test_has_labels_case_sensitivity(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        assert labels_upper
        tns.is_case_sensitive = True
        self.assertFalse(tns.has_taxa_labels(labels_upper))
        tns.is_case_sensitive = False
        self.assertTrue(tns.has_taxa_labels(labels_upper))

    def test_findall_multiple(self):
        tns = TaxonNamespace(self.str_labels)
        multilabels= ["_", "z"]
        for label in multilabels:
            tns.is_case_sensitive=True
            taxa = tns.findall(label=label)
            self.assertTrue(isinstance(taxa, collections.Iterable))
            self.assertEqual(len(taxa), len([s for s in self.str_labels if s == label]))
            for t in taxa:
                self.assertEqual(t.label, label)

    def test_findall_not_found(self):
        tns = TaxonNamespace(self.str_labels)
        tns.is_case_sensitive=True
        taxa = tns.findall(label="x")
        self.assertEqual(taxa, [])

    def test_get_taxon_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        for label in self.str_labels:
            t = tns.get_taxon(label)
            self.assertEqual(t.label, label)

    def test_get_nonexistant_taxon_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        check = ["u", "x", "y",]
        for label in check:
            assert label not in self.str_labels
            t = tns.get_taxon(check)
            self.assertIs(t, None)

    def test_case_insensitive_get_taxon_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        assert labels_upper
        # default: case insensitive
        for label in labels_upper:
            t = tns.get_taxon(label)
            self.assertIsNot(t, None)
            self.assertEqual(t.label.lower(), label.lower())
        # test: case sensitive
        tns.is_case_sensitive = True
        for label in labels_upper:
            t = tns.get_taxon(label)
            self.assertIs(t, None)

    def test_require_taxon_by_label_noadd(self):
        tns = TaxonNamespace(self.str_labels)
        for label in self.str_labels:
            t = tns.get_taxon(label)
            self.assertEqual(t.label, label)
        self.assertEqual(len(tns), len(self.str_labels))
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_require_taxon_by_label_add(self):
        tns = TaxonNamespace(self.str_labels)
        check = ["u", "x", "y",]
        for label in check:
            assert label not in self.str_labels
            t = tns.require_taxon(label)
            self.assertTrue(isinstance(t, Taxon))
            self.assertEqual(t.label, label)
        total = self.str_labels + check
        self.assertEqual(len(tns), len(total))
        self.validate_taxon_concepts(tns, total)

    def test_case_insensitive_require_taxon_by_label1(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        assert labels_upper
        for label in labels_upper:
            tns.is_case_sensitive = False
            t = tns.require_taxon(label)
            self.assertEqual(t.label.lower(), label.lower())
            self.assertEqual(len(tns), len(self.str_labels))
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_case_insensitive_require_taxon_by_label2(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        labels_upper = list(set(labels_upper))
        assert labels_upper
        for label in labels_upper:
            tns.is_case_sensitive = True
            t = tns.require_taxon(label)
            self.assertEqual(t.label, label)
        self.validate_taxon_concepts(tns, self.str_labels + labels_upper)

    def test_require_taxon_by_label_add_to_immutable(self):
        tns = TaxonNamespace(self.str_labels)
        tns.is_mutable = False
        check = ["u", "x", "y",]
        for label in check:
            assert label not in self.str_labels
            with self.assertRaises(TypeError):
                t = tns.require_taxon(label)

    def test_get_taxa_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        # label_set = set(self.str_labels)
        # taxa = tns.get_taxa(label_set)
        taxa = tns.get_taxa(self.str_labels + ["u", "x", "y"])
        self.assertEqual(len(taxa), len(self.str_labels))
        tx = [t.label for t in taxa]
        self.assertEqual(tx, self.str_labels)

    def test_get_nonexistant_taxa_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        check = ["u", "x", "y",]
        taxa = tns.get_taxa(check)
        self.assertEqual(len(taxa), 0)

    def test_case_insensitive_get_taxa_by_label(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        assert labels_upper
        # default: case-insensitive
        t2 = tns.get_taxa(labels_upper)
        self.assertEqual(len(t2), len(labels_upper))
        for t, label in zip(t2, labels_upper):
            self.assertEqual(t.label.lower(), label.lower())
        # test: case sensitive
        tns.is_case_sensitive = True
        t1 = tns.get_taxa(labels_upper)
        self.assertEqual(len(t1), 0)

    ### iteration ###

    def test_iter1(self):
        tns = TaxonNamespace(self.str_labels)
        for t1, label in zip(tns, self.str_labels):
            self.assertEqual(t1.label, label)

    def test_iter2(self):
        tns = TaxonNamespace(self.str_labels)
        for idx, t1 in enumerate(tns):
            self.assertEqual(t1.label, self.str_labels[idx])

    def test_reverse_iter(self):
        tns = TaxonNamespace(self.str_labels)
        r = self.str_labels[:]
        r.reverse()
        assert r != self.str_labels
        for idx, t1 in enumerate(reversed(tns)):
            self.assertEqual(t1.label, r[idx])

    ### sorting ###

    def test_sort(self):
        r = self.str_labels[:]
        r.sort()
        r.reverse()
        tns = TaxonNamespace(r)
        tns.sort()
        r2 = sorted(r)
        assert r != r2
        for idx, t1 in enumerate(tns):
            self.assertEqual(t1.label, r2[idx])

    def test_reverse(self):
        r = self.str_labels[:]
        r.sort()
        tns = TaxonNamespace(r)
        tns.reverse()
        r2 = r[:]
        r2.reverse()
        assert r != r2
        for idx, t1 in enumerate(tns):
            self.assertEqual(t1.label, r2[idx])

    def test_sorted(self):
        r = self.str_labels[:]
        r.sort()
        r.reverse()
        tns = TaxonNamespace(r)
        r2 = sorted(r)
        assert r != r2
        for idx, t1 in enumerate(sorted(tns)):
            self.assertEqual(t1.label, r2[idx])

    def test_reversed(self):
        r = self.str_labels[:]
        r.sort()
        tns = TaxonNamespace(r)
        r2 = r[:]
        r2.reverse()
        assert r != r2
        for idx, t1 in enumerate(reversed(tns)):
            self.assertEqual(t1.label, r2[idx])

    ### deletion ###

    def test_delete_by_index(self):
        for idx in range(len(self.taxa)):
            tns = TaxonNamespace(self.taxa)
            del tns[idx]
            for idx2, taxon in enumerate(self.taxa):
                if idx2 == idx:
                    self.assertNotIn(taxon, tns)
                else:
                    self.assertIn(taxon, tns)

    def test_remove_taxon(self):
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        for idx, taxon in enumerate(taxa):
            tns.remove_taxon(taxon)
            expected.remove(taxon)
            self.assertEqual(len(tns), len(expected))
            for idx2, taxon2 in enumerate(expected):
                if taxon2 in expected:
                    self.assertIn(taxon2, tns)
                elif taxon2 not in expected:
                    self.assertNotIn(taxon2, tns)

    def test_remove_taxon_error(self):
        tns = TaxonNamespace(self.str_labels)
        with self.assertRaises(ValueError):
            tns.remove_taxon(self.taxa[0])

    def test_remove_taxon_label(self):
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        for idx, label in enumerate(set(self.str_labels)):
            tns.remove_taxon_label(label)
            for t in taxa:
                if t.label == label and t in expected:
                    expected.remove(t)
            self.assertEqual(len(tns), len(expected))
            for t1, t2 in zip(tns, expected):
                self.assertIs(t1, t2)

    def test_remove_taxon_label_error(self):
        tns = TaxonNamespace(self.str_labels)
        key = "zzz"
        assert key not in self.str_labels
        with self.assertRaises(LookupError):
            tns.remove_taxon_label(key)

    def test_remove_taxon_label_case_insensitive(self):
        ucase_labels = [s.upper() for s in self.str_labels]
        assert ucase_labels
        assert ucase_labels != self.str_labels
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        for idx, label in enumerate(set(ucase_labels)):
            if label != label.lower():
                with self.assertRaises(LookupError):
                    tns.is_case_sensitive = True
                    tns.remove_taxon_label(label)
            tns.is_case_sensitive = False
            tns.remove_taxon_label(label)
            for t in taxa:
                if t.label.upper() == label.upper() and t in expected:
                    expected.remove(t)
            self.assertEqual(len(tns), len(expected))
            for t1, t2 in zip(tns, expected):
                self.assertIs(t1, t2)

    def test_discard_taxon_label(self):
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        for idx, label in enumerate(set(self.str_labels)):
            tns.discard_taxon_label(label)
            for t in taxa:
                if t.label == label and t in expected:
                    expected.remove(t)
            self.assertEqual(len(tns), len(expected))
            for t1, t2 in zip(tns, expected):
                self.assertIs(t1, t2)

    def test_discard_taxon_label_error(self):
        tns = TaxonNamespace(self.str_labels)
        key = "zzz"
        assert key not in self.str_labels
        try:
            tns.discard_taxon_label(key)
        except LookupError:
            self.fail()
        else:
            self.validate_taxon_concepts(tns, self.str_labels)

    def test_discard_taxon_label_case_insensitive(self):
        ucase_labels = [s.upper() for s in self.str_labels]
        assert ucase_labels
        assert ucase_labels != self.str_labels
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        # default: case-insensitive
        for idx, label in enumerate(set(ucase_labels)):
            tns.discard_taxon_label(label)
            for t in taxa:
                if t.label.upper() == label.upper() and t in expected:
                    expected.remove(t)
            self.assertEqual(len(tns), len(expected))
            for t1, t2 in zip(tns, expected):
                self.assertIs(t1, t2)

    def test_discard_taxon_label_case_sensitive(self):
        ucase_labels = [s.upper() for s in self.str_labels]
        assert ucase_labels
        assert ucase_labels != self.str_labels
        taxa = [Taxon(s) for s in self.str_labels]
        tns = TaxonNamespace(taxa)
        expected = taxa[:]
        # test: case sensitive
        tns.is_case_sensitive = True
        for idx, label in enumerate(set(ucase_labels)):
            if label != label.lower():
                x1 = len(tns)
                try:
                    tns.discard_taxon_label(label)
                except LookupError:
                    self.fail()
                else:
                    self.assertEqual(len(tns), x1)

    def test_clear(self):
        tns = TaxonNamespace(self.str_labels)
        self.assertEqual(len(tns), len(self.str_labels))
        tns.clear()
        self.assertEqual(len(tns), 0)
        x = []
        for t in tns:
            x.append(t)
        self.assertEqual(len(x), 0)

class TaxonNamespaceIdentity(unittest.TestCase):

    def setUp(self):
        self.str_labels = ["a", "a", "b", "c", "d", "e", "_", "_", "_", "z", "z", "z"]
        self.taxa = [ Taxon(label) for label in self.str_labels ]
        self.tns1 = TaxonNamespace(self.taxa)
        self.tns2 = TaxonNamespace(self.taxa)

    # def test_separate_but_equal(self):
    #     self.assertIsNot(self.tns1, self.tns2)
    #     self.assertEqual(self.tns1, self.tns2)

    # def test_different_labels(self):
    #     self.assertIsNot(self.tns1, self.tns2)
    #     self.assertEqual(self.tns1, self.tns2)
    #     self.tns1.label = "hello"
    #     self.tns2.label = "goodbye"
    #     self.assertNotEqual(self.tns1, self.tns2)
    #     self.tns1.label = self.tns2.label
    #     self.assertIsNot(self.tns1, self.tns2)
    #     self.assertEqual(self.tns1, self.tns2)

    # def test_same_annotations(self):
    #     self.assertIsNot(self.tns1, self.tns2)
    #     self.assertEqual(self.tns1, self.tns2)
    #     self.tns1.annotations.add_new("hello", 0)
    #     self.tns2.annotations.add_new("hello", 0)
    #     # not equal because each |AnnotationSet| has a ``target``
    #     # attribute that holds reference to the object being annotated. As
    #     # these the target objects are necessarily different (even if they
    #     # evaluate being equal), the |AnnotationSet| objects are not
    #     # considered equal, and thus the target objects that have the
    #     # AnnotationSets are not equal.
    #     self.assertNotEqual(self.tns1, self.tns2)

    # def test_different_annotations1(self):
    #     self.assertIsNot(self.tns1, self.tns2)
    #     self.assertEqual(self.tns1, self.tns2)
    #     self.tns1.annotations.add_new("hello", 0)
    #     self.assertNotEqual(self.tns1, self.tns2)

    def test_hash_dict_membership(self):
        k = {}
        k[self.tns1] = 1
        k[self.tns2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.tns1], 1)
        self.assertEqual(k[self.tns2], 2)
        self.assertIn(self.tns1, k)
        self.assertIn(self.tns2, k)
        del k[self.tns1]
        self.assertNotIn(self.tns1, k)
        self.assertIn(self.tns2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.tns1: 1}
        k2 = {self.tns2: 1}
        self.assertIn(self.tns1, k1)
        self.assertIn(self.tns2, k2)
        self.assertNotIn(self.tns2, k1)
        self.assertNotIn(self.tns1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.tns1)
        k.add(self.tns2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.tns1, k)
        self.assertIn(self.tns2, k)
        k.discard(self.tns1)
        self.assertNotIn(self.tns1, k)
        self.assertIn(self.tns2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.tns1: 1}
        k2 = {self.tns2: 1}
        self.assertIn(self.tns1, k1)
        self.assertIn(self.tns2, k2)
        self.assertNotIn(self.tns2, k1)
        self.assertNotIn(self.tns1, k2)

# note that compare_and_validate.Comparator must be listed first,
# otherwise setUp will not be called
class TaxonNamespaceCloning(compare_and_validate.Comparator, unittest.TestCase):

    def setUp(self):
        self.str_labels = ["a", "a", "b", "c", "d", "e", "_", "_", "_", "z", "z", "z"]
        self.taxa = [ Taxon(label) for label in self.str_labels ]
        self.tns1 = TaxonNamespace(self.taxa, label="T1")

    def test_taxon_namespace_scoped_copy(self):
        for tns2 in (self.tns1.clone(1),
                self.tns1.taxon_namespace_scoped_copy()):
            self.assertIs(tns2, self.tns1)

    def test_construct_from_another(self):
        tns2 = TaxonNamespace(self.tns1)
        for tns2 in (TaxonNamespace(self.tns1),
                self.tns1.clone(0),
                copy.copy(self.tns1)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2.label, self.tns1.label)
            self.assertEqual(tns2._taxa, self.tns1._taxa)
            for t1, t2 in zip(self.tns1, tns2):
                self.assertIs(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_construct_from_another_different_label(self):
        tns2 = TaxonNamespace(self.tns1, label="T2")
        self.assertIsNot(tns2, self.tns1)
        self.assertNotEqual(tns2.label, self.tns1.label)
        self.assertEqual(self.tns1.label, "T1")
        self.assertEqual(tns2.label, "T2")
        self.assertEqual(tns2._taxa, self.tns1._taxa)
        for t1, t2 in zip(self.tns1, tns2):
            self.assertIs(t1, t2)
        self.compare_distinct_annotables(tns2, self.tns1)

    def test_construct_from_changed_label(self):
        for tns2 in (TaxonNamespace(self.tns1),
                self.tns1.clone(0),
                copy.copy(self.tns1)):
            tns2.label = "T2"
            self.assertNotEqual(tns2.label, self.tns1.label)
            self.assertEqual(self.tns1.label, "T1")
            self.assertEqual(tns2.label, "T2")
            self.assertEqual(tns2._taxa, self.tns1._taxa)
            for t1, t2 in zip(self.tns1, tns2):
                self.assertIs(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_construct_from_another_with_simple_annotations(self):
        self.tns1.annotations.add_new("A", 1)
        self.tns1.annotations.add_new("B", 2)
        self.tns1.annotations.add_new("C", 3)
        for tns2 in (TaxonNamespace(self.tns1),
                self.tns1.clone(0),
                copy.copy(self.tns1)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2._taxa, self.tns1._taxa)
            for t1, t2 in zip(tns2, self.tns1):
                self.assertIs(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_construct_from_another_with_complex_annotations(self):
        self.tns1.annotations.add_new("a", 0)
        b = self.tns1.annotations.add_new("b", (self.tns1, "label"), is_attribute=True)
        b.annotations.add_new("c", 3)
        self.tns1.annotations.add_new("A", 1)
        self.tns1.annotations.add_new("B", 2)
        self.tns1.annotations.add_new("C", 3)
        for tns2 in (TaxonNamespace(self.tns1), self.tns1.clone(0), copy.copy(self.tns1)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2._taxa, self.tns1._taxa)
            for t1, t2 in zip(tns2, self.tns1):
                self.assertIs(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_deepcopy_from_another(self):
        for tns2 in (copy.deepcopy(self.tns1),
                self.tns1.clone(2)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2.label, self.tns1.label)
            self.assertEqual(len(tns2), len(self.tns1))
            for t1, t2 in zip(self.tns1, tns2):
                self.assertIsNot(t1, t2)
                self.assertEqual(t1.label, t2.label)
                self.compare_distinct_annotables(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_deepcopy_from_another_with_simple_annotations(self):
        self.tns1.annotations.add_new("a", 0)
        self.tns1.annotations.add_new("b", 1)
        self.tns1.annotations.add_new("c", 3)
        for tns2 in (copy.deepcopy(self.tns1),
                self.tns1.clone(2)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2.label, self.tns1.label)
            self.assertEqual(len(tns2), len(self.tns1))
            for t1, t2 in zip(self.tns1, tns2):
                self.assertIsNot(t1, t2)
                self.assertEqual(t1.label, t2.label)
                self.compare_distinct_annotables(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)

    def test_deepcopy_from_another_with_complex_annotations(self):
        self.tns1.annotations.add_new("a", 0)
        b = self.tns1.annotations.add_new("b", (self.tns1, "label"), is_attribute=True)
        b.annotations.add_new("c", 3)
        for tns2 in (copy.deepcopy(self.tns1),
                self.tns1.clone(2)):
            self.assertIsNot(tns2, self.tns1)
            self.assertEqual(tns2.label, self.tns1.label)
            self.assertEqual(len(tns2), len(self.tns1))
            for t1, t2 in zip(self.tns1, tns2):
                self.assertIsNot(t1, t2)
                self.assertEqual(t1.label, t2.label)
                self.compare_distinct_annotables(t1, t2)
            self.compare_distinct_annotables(tns2, self.tns1)
            self.tns1.label = "x"
            tns2.label = "y"
            self.assertEqual(self.tns1.annotations[1].value, "x")
            self.assertEqual(tns2.annotations[1].value, "y")
            self.tns1.label = "T1"

if __name__ == "__main__":
    unittest.main()

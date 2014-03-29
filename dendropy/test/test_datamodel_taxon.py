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
Tests creation, reading, update, deletion of Taxon and TaxonNamespace objects.
"""

import collections
import unittest
from dendropy import Taxon, TaxonNamespace

class TaxonNamespaceLinkedTest(unittest.TestCase):

    pass

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
        self.assertEqual(tns1, tns2)

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
            self.assertFalse(tns.has_taxon_label(label))
            self.assertTrue(tns.has_taxon_label(label, case_insensitive=True))

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
        self.assertFalse(tns.has_taxa_labels(labels_upper))
        self.assertTrue(tns.has_taxa_labels(labels_upper, case_insensitive=True))

    def test_findall_multiple(self):
        tns = TaxonNamespace(self.str_labels)
        multilabels= ["_", "z"]
        for label in multilabels:
            taxa = tns.findall(label=label, case_insensitive=False)
            self.assertTrue(isinstance(taxa, collections.Iterable))
            self.assertEqual(len(taxa), len([s for s in self.str_labels if s == label]))
            for t in taxa:
                self.assertEqual(t.label, label)

    def test_findall_not_found(self):
        tns = TaxonNamespace(self.str_labels)
        taxa = tns.findall(label="x", case_insensitive=False)
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
        for label in labels_upper:
            t = tns.get_taxon(label)
            self.assertIs(t, None)
            t = tns.get_taxon(label, case_insensitive=True)
            self.assertIsNot(t, None)
            self.assertEqual(t.label.lower(), label.lower())

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
            t = tns.require_taxon(label, case_insensitive=True)
            self.assertEqual(t.label.lower(), label.lower())
            self.assertEqual(len(tns), len(self.str_labels))
        self.validate_taxon_concepts(tns, self.str_labels)

    def test_case_insensitive_require_taxon_by_label2(self):
        tns = TaxonNamespace(self.str_labels)
        labels_upper = [label.upper() for label in self.str_labels if label.upper() != label]
        labels_upper = list(set(labels_upper))
        assert labels_upper
        for label in labels_upper:
            t = tns.require_taxon(label, case_insensitive=False)
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
        t1 = tns.get_taxa(labels_upper)
        self.assertEqual(len(t1), 0)
        t2 = tns.get_taxa(labels_upper, case_insensitive=True)
        self.assertEqual(len(t2), len(labels_upper))
        for t, label in zip(t2, labels_upper):
            self.assertEqual(t.label.lower(), label.lower())

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
                    tns.remove_taxon_label(label)
            tns.remove_taxon_label(label, case_insensitive=True)
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
        for idx, label in enumerate(set(ucase_labels)):
            if label != label.lower():
                x1 = len(tns)
                try:
                    tns.discard_taxon_label(label)
                except LookupError:
                    self.fail()
                else:
                    self.assertEqual(len(tns), x1)
            tns.discard_taxon_label(label, case_insensitive=True)
            for t in taxa:
                if t.label.upper() == label.upper() and t in expected:
                    expected.remove(t)
            self.assertEqual(len(tns), len(expected))
            for t1, t2 in zip(tns, expected):
                self.assertIs(t1, t2)

    def test_clear(self):
        tns = TaxonNamespace(self.str_labels)
        self.assertEqual(len(tns), len(self.str_labels))
        tns.clear()
        self.assertEqual(len(tns), 0)
        x = []
        for t in tns:
            x.append(t)
        self.assertEqual(len(x), 0)

    # def test_drop_taxon_singlekeyed_singledropped_casesensitive_noerror(self):
    #     upper_labels = [label.upper() for label in self.str_labels if label.upper() != label]
    #     full_list = self.str_labels + upper_labels
    #     tns = TaxonNamespace(full_list)
    #     for label in full_list:
    #         tns.drop_taxon(key=label,
    #                 multiple=False,
    #                 case_insensitive=False,
    #                 error_if_not_found=False)
    #         z = full_list[:]
    #         z.remove(label)
    #         self.validate_taxon_concepts(tns, z)


#     def test_delete_by_label(self):
#         for idx in range(len(self.str_labels)):
#             tns = TaxonNamespace(self.str_labels)
#             expected_labels = [label in self.str_labels if label != self.str_labels[idx]]
#             del tns[idx]
#             for idx2, taxon in enumerate(self.taxa):
#                 if idx2 == idx:
#                     self.assertNotIn(taxon, tns)
#                 else:
#                     self.assertIn(taxon, tns)



# class TaxaTest(datatest.AnnotatedDataObjectVerificationTestCase):

#     def setUp(self):
#         self.labels = []
#         for idx in xrange(10):
#             self.labels.append("T%d" % (idx+1))
#         self.taxon_set = dendropy.TaxonSet()
#         for label in self.labels:
#             self.taxon_set.new_taxon(label=label)

#     def testLabelsAsKeys(self):
#         for t in self.taxon_set:
#             self.assertIs(t, self.taxon_set[t.label])

#     def testPositiveIndexing(self):
#         for i, t in enumerate(self.taxon_set):
#             self.assertIs(t, self.taxon_set[i])

#     def testNegativeIndexing(self):
#         for i, t in enumerate(self.taxon_set):
#             self.assertIs(t, self.taxon_set[i])

#     def testRaisesKeyError(self):
#         self.assertRaises(KeyError, self.taxon_set.__getitem__, 'foo')
#         self.assertRaises(KeyError, self.taxon_set.__getitem__, 'T')

#     def testRaisesIndexError(self):
#         self.assertRaises(IndexError, self.taxon_set.__getitem__, 1000)
#         self.assertRaises(IndexError, self.taxon_set.__getitem__, -1000)

#     def testCompositionFromStrings(self):
#         ts = dendropy.TaxonSet(self.labels)
#         self.assertDistinctButEqual(ts, self.taxon_set)

#     def testCompositionFromTaxa(self):
#         ts = dendropy.TaxonSet(self.taxon_set)
#         self.assertDistinctButEqual(ts, self.taxon_set, distinct_taxon_objects=False)

#     def testTaxaQuerying(self):
#         ts = dendropy.TaxonSet(self.labels)
#         self.assertTrue(ts.has_taxa(labels=self.labels))
#         self.assertTrue(ts.has_taxa(taxa=ts))
#         self.assertFalse(ts.has_taxa(labels=self.labels+["k"]))
#         k = ts.new_taxon(label="k")
#         self.assertTrue(ts.has_taxa(taxa=[k]))
#         self.assertTrue(ts.has_taxon(label="k"))
#         self.assertTrue(ts.has_taxa(labels=self.labels+["k"]))
#         j = dendropy.Taxon(label="j")
#         ts.add_taxon(j)
#         self.assertTrue(ts.has_taxa(taxa=[j]))
#         self.assertTrue(ts.has_taxon(label="j"))
#         self.assertTrue(ts.has_taxa(labels=self.labels+["j"]))
#         self.assertFalse(ts.has_taxon(taxon=dendropy.Taxon()))
#         for label in self.labels:
#             self.assertTrue(ts.has_taxon(label=label))

#     def testLockedVsUnlocked(self):
#         self.taxon_set.lock()
#         self.assertEquals(len(self.taxon_set), 10)
#         for idx, t in enumerate(self.taxon_set):
#             self.assertEquals(t.label, self.labels[idx])
#         self.assertRaises(KeyError, self.taxon_set.new_taxon, label="A1")
#         self.assertRaises(KeyError, self.taxon_set.require_taxon, label="A1", oid=None)
#         self.taxon_set.unlock()
#         x1 = self.taxon_set.new_taxon(label="X1")
#         self.assertIs(x1, self.taxon_set.get_taxon(label="X1"))
#         self.assertIs(self.taxon_set.get_taxon(label="X2"), None)
#         self.taxon_set.require_taxon(label="X3")
#         self.assertEquals(len(self.taxon_set), 12)

#     def testTaxonQuerying(self):
#         ts = dendropy.TaxonSet(self.labels)
#         self.assertIs(ts.get_taxon(label="Q"), None)
#         self.assertIs(ts.get_taxon(label="T1"), ts[0])

# class TaxonSetPartitionTest(datatest.AnnotatedDataObjectVerificationTestCase):

#     def setUp(self):
#         self.taxon_set = dendropy.TaxonSet([
#                 'a1', 'a2', 'a3', 'a4',
#                 'b1', 'b2', 'b3', 'b4',
#                 'c1', 'c2', 'c2', 'c3',
#                 'd1', 'a5', 'a6', 'd2',
#                 'd3'])
#         self.membership_func = lambda x: x.label[0]
#         self.membership_dict = {}
#         for t in self.taxon_set:
#             self.membership_dict[t] = t.label[0]
#             t.subset_id = t.label[0]
#         self.membership_lists = [
#             [self.taxon_set[0], self.taxon_set[1], self.taxon_set[2], self.taxon_set[3],
#              self.taxon_set[13], self.taxon_set[14]],
#             [self.taxon_set[4], self.taxon_set[5], self.taxon_set[6], self.taxon_set[7]],
#             [self.taxon_set[8], self.taxon_set[9], self.taxon_set[10], self.taxon_set[11]],
#             [self.taxon_set[12], self.taxon_set[15], self.taxon_set[16]]]
#         self.label_map = ['a', 'b', 'c', 'd']
#         self.expected_sets = set([dendropy.TaxonSet(s, label=self.label_map[i]) \
#                 for i, s in enumerate(self.membership_lists)])
#         self.expected_dict = {}
#         for s in self.expected_sets:
#             self.expected_dict[self.membership_dict[s[0]]] = s

#     def verify_subsets(self, subsets, use_label_indices=False):
#         for s in subsets:
#             if use_label_indices:
#                 key = self.label_map[s.label]
#             else:
#                 key = s.label
#             self.assertDistinctButEqual(
#                     self.expected_dict[key],
#                     s,
#                     distinct_taxa=True,
#                     distinct_taxon_objects=False)

#     def testFromMembershipFunc(self):
#         tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_func=self.membership_func)
#         self.verify_subsets(tsp.subsets())

#     def testFromMembershipAttr(self):
#         tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_attr_name='subset_id')
#         self.verify_subsets(tsp.subsets())

#     def testFromMembershipDict(self):
#         tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_dict=self.membership_dict)
#         self.verify_subsets(tsp.subsets())

#     def testFromMembershipLists(self):
#         tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_lists=self.membership_lists)
#         self.verify_subsets(tsp.subsets(), use_label_indices=True)

# class TaxonSetMappingTest(datatest.AnnotatedDataObjectVerificationTestCase):

#     def setUp(self):
#         self.domain_taxa = dendropy.TaxonSet([
#                 'a1', 'a2', 'a3', 'a4',
#                 'b1', 'b2', 'b3', 'b4',
#                 'c1', 'c2', 'c2', 'c3',
#                 'd1', 'a5', 'a6', 'd2',
#                 'd3'])
#         self.range_taxa = dendropy.TaxonSet([
#             'A', 'B', 'C', 'D',])
#         self.domain_taxa.lock()
#         self.range_taxa.lock()
#         self.mapping_func = lambda x: self.range_taxa.require_taxon(label=x.label[0].upper())
#         self.mapping_dict = {}
#         for t in self.domain_taxa:
#             self.mapping_dict[t] = self.mapping_func(t)
#             t.containing_taxa = self.mapping_dict[t]
#         self.expected_forward_label_map = {
#                 'a1' : 'A',
#                 'a2' : 'A',
#                 'a3' : 'A',
#                 'a4' : 'A',
#                 'a5' : 'A',
#                 'a6' : 'A',
#                 'b1' : 'B',
#                 'b2' : 'B',
#                 'b3' : 'B',
#                 'b4' : 'B',
#                 'c1' : 'C',
#                 'c2' : 'C',
#                 'c3' : 'C',
#                 'c4' : 'C',
#                 'd1' : 'D',
#                 'd2' : 'D',
#                 'd3' : 'D',}
#         self.expected_backward_label_map = {
#                 'A' : set(['a1', 'a2', 'a3', 'a4', 'a5', 'a6']),
#                 'B' : set(['b1', 'b2', 'b3', 'b4',]),
#                 'C' : set(['c1', 'c2', 'c2', 'c3',]),
#                 'D' : set(['d1', 'd2', 'd3'])
#                 }

#     def verifyMapping(self, tsm):
#         for t in self.domain_taxa:
#             self.assertEqual(tsm.forward[t].label, self.expected_forward_label_map[t.label])
#         for t in self.range_taxa:
#             self.assertEqual(set([i.label for i in tsm.reverse[t]]), self.expected_backward_label_map[t.label])

#     def testFromFunc(self):
#         tsm = dendropy.TaxonSetMapping(mapping_func=self.mapping_func, domain_taxon_set=self.domain_taxa)
#         self.verifyMapping(tsm)

#     def testFromAttr(self):
#         tsm = dendropy.TaxonSetMapping(mapping_attr_name='containing_taxa', domain_taxon_set=self.domain_taxa)
#         self.verifyMapping(tsm)

#     def testFromDict(self):
#         tsm = dendropy.TaxonSetMapping(mapping_dict=self.mapping_dict)
#         self.verifyMapping(tsm)

# class FullCopyTaxaTestCase(datatest.AnnotatedDataObjectVerificationTestCase):

#     def testFullCopyTaxonSet(self):
#         src_path = pathmap.char_source_path("crotaphytus_bicinctores.cytb.aligned.nexml")
#         d = dendropy.DataSet.get_from_path(src_path, "nexml")
#         taxon_set1 = d.taxon_sets[0]
#         taxon_set2 = taxon_set1.fullcopy()
#         self.assertTrue(taxon_set1 is not taxon_set2)
#         self.assertDistinctButEqualTaxonSet(taxon_set1, taxon_set2, distinct_taxon_objects=True)
#         self.assertEqual(len(taxon_set1), len(taxon_set2))
#         for idx, taxon1 in enumerate(taxon_set1):
#             taxon2 = taxon_set2[idx]
#             self.assertTrue(taxon1 is not taxon2)
#             self.assertEqual(taxon1.label, taxon2.label)
#             for idx, a1 in enumerate(taxon1.annotations):
#                 a2 = taxon2.annotations[idx]
#                 self.assertTrue(a1 is not a2)
#                 self.assertEqual(a1.name, a2.name)
#                 self.assertEqual(a1.value, a2.value)
#         self.assertEqual(len(taxon_set1.annotations), len(taxon_set2.annotations))
#         for idx, a1 in enumerate(taxon_set1.annotations):
#             a2 = taxon_set2.annotations[idx]
#             self.assertTrue(a1 is not a2)
#             self.assertEqual(a1.name, a2.name)
#             self.assertEqual(a1.value, a2.value)

if __name__ == "__main__":
    unittest.main()

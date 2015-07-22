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
Tests for NEWICK taxon handling.
"""

import sys
import os
import unittest
import dendropy
from dendropy.utility import error
from dendropy.test.support import dendropytest
from dendropy.dataio import nexusreader
from dendropy.dataio import nexusprocessing

class TaxonSymbolMappingTest(unittest.TestCase):

    def test_standard_lookup_and_create(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        for idx, label in enumerate(labels):
            self.assertEqual(len(tns), idx)
            t1 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertEqual(t1.label, label)

            t2 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t2)
            self.assertEqual(t2.label, label)

            t3 = tsm.require_taxon_for_symbol(str(idx+1))
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t3)
            self.assertEqual(t3.label, label)

            assert label.upper() != label
            t4 = tsm.require_taxon_for_symbol(label.upper())
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t4)
            self.assertEqual(t4.label, label)

    def test_no_number_lookup_and_create(self):
        # looking up a number symbol should result in new taxon creation
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns,
                enable_lookup_by_taxon_number=False)
        for idx, label in enumerate(labels):
            self.assertEqual(len(tns), idx)
            t1 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertEqual(t1.label, label)

            t2 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t2)
            self.assertEqual(t2.label, label)

            t3 = tsm.lookup_taxon_symbol(str(idx+1), create_taxon_if_not_found=False)
            self.assertIs(t3, None)
            self.assertEqual(len(tns), idx+1)

    def test_no_number_lookup_and_create2(self):
        # looking up a number symbol should result in new taxon creation
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns,
                enable_lookup_by_taxon_number=False)
        taxa = []
        for label_idx, label in enumerate(labels):
            t = dendropy.Taxon(label)
            tsm.add_taxon(t)
            taxa.append(t)
        self.assertEqual(len(tns), len(labels))
        for label_idx, label in enumerate(labels):
            t1 = tsm.require_taxon_for_symbol(label_idx+1)
            self.assertNotIn(t1, taxa)
            self.assertEqual(t1.label, str(label_idx+1))
            self.assertEqual(len(tns), len(labels)+label_idx+1)

    def test_new_taxon(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        for label_idx, label in enumerate(labels):
            t = tsm.new_taxon(label)
            self.assertEqual(len(tns), label_idx+1)
            self.assertEqual(t.label, label)
            self.assertIs(tsm.require_taxon_for_symbol(label), t)
            self.assertEqual(len(tns), label_idx+1)
            self.assertIs(tsm.require_taxon_for_symbol(str(label_idx+1)), t)
            self.assertEqual(len(tns), label_idx+1)
        self.assertEqual(len(tns), len(labels))

    def test_add_taxon(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        for label_idx, label in enumerate(labels):
            t = dendropy.Taxon(label)
            tsm.add_taxon(t)
            self.assertEqual(len(tns), label_idx+1)
            self.assertEqual(t.label, label)
            self.assertIs(tsm.require_taxon_for_symbol(label), t)
            self.assertEqual(len(tns), label_idx+1)
            self.assertIs(tsm.require_taxon_for_symbol(str(label_idx+1)), t)
            self.assertEqual(len(tns), label_idx+1)
        self.assertEqual(len(tns), len(labels))

    def test_simple_token_lookup(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        translate = {}
        t_labels = {}
        for label_idx, label in enumerate(labels):
            t = dendropy.Taxon(label)
            t_labels[t] = t.label
            tsm.add_taxon(t)
            token = label_idx + 1
            translate[token] = t
            tsm.add_translate_token(token, t)
        self.assertEqual(len(tns), len(labels))
        for token in translate:
            t1 = translate[token]
            t2 = tsm.require_taxon_for_symbol(token)
            self.assertIs(t1, t2)
            self.assertEqual(t2.label, t_labels[t1])
        self.assertEqual(len(tns), len(labels))

    def test_tricky_token_lookup(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        translate = {}
        t_labels = {}
        for label_idx, label in enumerate(labels):
            t = dendropy.Taxon(label)
            t_labels[t] = t.label
            tsm.add_taxon(t)
            token = str(len(labels) - label_idx)
            translate[token] = t
            tsm.add_translate_token(token, t)
        self.assertEqual(len(tns), len(labels))
        for token in translate:
            t1 = translate[token]
            t2 = tsm.require_taxon_for_symbol(token)
            self.assertIs(t1, t2)
            self.assertEqual(t2.label, t_labels[t1])
        self.assertEqual(len(tns), len(labels))

    def test_mixed_token_lookup(self):
        labels = ["t{}".format(i) for i in range(1, 101)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        translate = {}
        t_labels = {}
        labels_t = {}
        for label_idx, label in enumerate(labels):
            t = dendropy.Taxon(label)
            t_labels[t] = t.label
            labels_t[t.label] = t
            tsm.add_taxon(t)
            if label_idx % 2 == 0:
                token = str(label_idx+1)
                translate[token] = t
                tsm.add_translate_token(token, t)
        self.assertEqual(len(tns), len(labels))
        for label_idx, label in enumerate(labels):
            token = label_idx + 1
            t1 = tsm.require_taxon_for_symbol(token)
            self.assertEqual(len(tns), len(labels))
            self.assertEqual(t1.label, label)
            self.assertIs(t1, labels_t[label])
            if token in translate:
                self.assertIs(t1, translate[token])
        self.assertEqual(len(tns), len(labels))

    def test_taxon_namespace_locking(self):
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        self.assertFalse(tns.is_mutable)
        del tsm
        self.assertTrue(tns.is_mutable)

class NexusTaxaCaseInsensitivityTest(unittest.TestCase):

    def setUp(self):
        self.data_str = """\
            #NEXUS

            BEGIN TAXA;
                DIMENSIONS NTAX=5;
                TAXLABELS AAA BBB CCC DDD EEE;
            END;

            BEGIN CHARACTERS;
                DIMENSIONS  NCHAR=8;
                FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. INTERLEAVE;
                MATRIX
                    AAA ACGT
                    BBB ACGT
                    CCC ACGT
                    DDD ACGT
                    EEE ACGT

                    aaa ACGT
                    bbb ACGT
                    ccc ACGT
                    ddd ACGT
                    eee ACGT
                ;
            END;
            """

    def testCaseInsensitiveChars(self):
        d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus', case_sensitive_taxon_labels=False)
        expected = ["AAA", "BBB", "CCC", "DDD", "EEE"]
        observed = [t.label.upper() for t in d.taxon_namespace]
        for i, x in enumerate(expected):
            self.assertTrue(x in observed)
        for i, x in enumerate(observed):
            self.assertTrue(x in expected)
        self.assertEqual(len(d.taxon_namespace), 5)

    def testCaseSensitiveChars(self):
        #d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus', case_sensitive_taxon_labels=False)
        self.assertRaises(error.DataParseError,
                dendropy.DnaCharacterMatrix.get_from_string,
                self.data_str,
                'nexus',
                case_sensitive_taxon_labels=True)

    def testDefaultCaseSensitivityChars(self):
        d = dendropy.DnaCharacterMatrix.get_from_string(self.data_str, 'nexus')
        expected = ["AAA", "BBB", "CCC", "DDD", "EEE"]
        observed = [t.label.upper() for t in d.taxon_namespace]
        for i, x in enumerate(expected):
            self.assertTrue(x in observed)
        for i, x in enumerate(observed):
            self.assertTrue(x in expected)
        self.assertEqual(len(d.taxon_namespace), 5)

class NexusTooManyTaxaTest(
        dendropytest.ExtendedTestCase):

    def testTooManyTaxaNonInterleaved(self):
        data_str = """\
        #NEXUS
        BEGIN TAXA;
            DIMENSIONS NTAX=2;
            TAXLABELS AAA BBB CCC DDD EEE;
        END;
        """
        self.assertRaises(nexusreader.NexusReader.TooManyTaxaError,
                dendropy.DnaCharacterMatrix.get_from_string,
                data_str,
                'nexus')

if __name__ == "__main__":
    unittest.main()

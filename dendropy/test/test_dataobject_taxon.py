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
Tests creation, reading, update, deletion of Taxon and TaxonSet objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import datatest
from dendropy.test.support.extendedtest import ExtendedTestCase
import dendropy

class TaxaTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.labels = []
        for idx in xrange(10):
            self.labels.append("T%d" % (idx+1))
        self.taxon_set = dendropy.TaxonSet()
        for label in self.labels:
            self.taxon_set.new_taxon(label=label)

    def testLabelsAsKeys(self):
        for t in self.taxon_set:
            self.assertIs(t, self.taxon_set[t.label])

    def testPositiveIndexing(self):
        for i, t in enumerate(self.taxon_set):
            self.assertIs(t, self.taxon_set[i])

    def testNegativeIndexing(self):
        for i, t in enumerate(self.taxon_set):
            self.assertIs(t, self.taxon_set[i])

    def testRaisesKeyError(self):
        self.assertRaises(KeyError, self.taxon_set.__getitem__, 'foo')
        self.assertRaises(KeyError, self.taxon_set.__getitem__, 'T')

    def testRaisesIndexError(self):
        self.assertRaises(IndexError, self.taxon_set.__getitem__, 1000)
        self.assertRaises(IndexError, self.taxon_set.__getitem__, -1000)

    def testCompositionFromStrings(self):
        ts = dendropy.TaxonSet(self.labels)
        self.assertDistinctButEqual(ts, self.taxon_set)

    def testCompositionFromTaxa(self):
        ts = dendropy.TaxonSet(self.taxon_set)
        self.assertDistinctButEqual(ts, self.taxon_set, distinct_taxon_objects=False)

    def testTaxaQuerying(self):
        ts = dendropy.TaxonSet(self.labels)
        self.assertTrue(ts.has_taxa(labels=self.labels))
        self.assertTrue(ts.has_taxa(taxa=ts))
        self.assertFalse(ts.has_taxa(labels=self.labels+["k"]))
        k = ts.new_taxon(label="k")
        self.assertTrue(ts.has_taxa(taxa=[k]))
        self.assertTrue(ts.has_taxon(label="k"))
        self.assertTrue(ts.has_taxa(labels=self.labels+["k"]))
        j = dendropy.Taxon(label="j")
        ts.add_taxon(j)
        self.assertTrue(ts.has_taxa(taxa=[j]))
        self.assertTrue(ts.has_taxon(label="j"))
        self.assertTrue(ts.has_taxa(labels=self.labels+["j"]))
        self.assertFalse(ts.has_taxon(taxon=dendropy.Taxon()))
        for label in self.labels:
            self.assertTrue(ts.has_taxon(label=label))

    def testLockedVsUnlocked(self):
        self.taxon_set.lock()
        self.assertEquals(len(self.taxon_set), 10)
        for idx, t in enumerate(self.taxon_set):
            self.assertEquals(t.label, self.labels[idx])
        self.assertRaises(KeyError, self.taxon_set.new_taxon, label="A1")
        self.assertRaises(KeyError, self.taxon_set.require_taxon, label="A1", oid=None)
        self.taxon_set.unlock()
        x1 = self.taxon_set.new_taxon(label="X1")
        self.assertIs(x1, self.taxon_set.get_taxon(label="X1"))
        self.assertIs(self.taxon_set.get_taxon(label="X2"), None)
        self.taxon_set.require_taxon(label="X3")
        self.assertEquals(len(self.taxon_set), 12)

    def testTaxonQuerying(self):
        ts = dendropy.TaxonSet(self.labels)
        self.assertIs(ts.get_taxon(label="Q"), None)
        self.assertIs(ts.get_taxon(label="T1"), ts[0])

class TaxonSetPartitionTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.taxon_set = dendropy.TaxonSet([
                'a1', 'a2', 'a3', 'a4',
                'b1', 'b2', 'b3', 'b4',
                'c1', 'c2', 'c2', 'c3',
                'd1', 'a5', 'a6', 'd2',
                'd3'])
        self.membership_func = lambda x: x.label[0]
        self.membership_dict = {}
        for t in self.taxon_set:
            self.membership_dict[t] = t.label[0]
            t.subset_id = t.label[0]
        self.membership_lists = [
            [self.taxon_set[0], self.taxon_set[1], self.taxon_set[2], self.taxon_set[3],
             self.taxon_set[13], self.taxon_set[14]],
            [self.taxon_set[4], self.taxon_set[5], self.taxon_set[6], self.taxon_set[7]],
            [self.taxon_set[8], self.taxon_set[9], self.taxon_set[10], self.taxon_set[11]],
            [self.taxon_set[12], self.taxon_set[15], self.taxon_set[16]]]
        self.label_map = ['a', 'b', 'c', 'd']
        self.expected_sets = set([dendropy.TaxonSet(s, label=self.label_map[i]) \
                for i, s in enumerate(self.membership_lists)])
        self.expected_dict = {}
        for s in self.expected_sets:
            self.expected_dict[self.membership_dict[s[0]]] = s

    def verify_subsets(self, subsets, use_label_indices=False):
        for s in subsets:
            if use_label_indices:
                key = self.label_map[s.label]
            else:
                key = s.label
            self.assertDistinctButEqual(
                    self.expected_dict[key],
                    s,
                    distinct_taxa=True,
                    distinct_taxon_objects=False)

    def testFromMembershipFunc(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_func=self.membership_func)
        self.verify_subsets(tsp.subsets())

    def testFromMembershipAttr(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_attr_name='subset_id')
        self.verify_subsets(tsp.subsets())

    def testFromMembershipDict(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_dict=self.membership_dict)
        self.verify_subsets(tsp.subsets())

    def testFromMembershipLists(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_lists=self.membership_lists)
        self.verify_subsets(tsp.subsets(), use_label_indices=True)

class TaxonSetMappingTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.domain_taxa = dendropy.TaxonSet([
                'a1', 'a2', 'a3', 'a4',
                'b1', 'b2', 'b3', 'b4',
                'c1', 'c2', 'c2', 'c3',
                'd1', 'a5', 'a6', 'd2',
                'd3'])
        self.range_taxa = dendropy.TaxonSet([
            'A', 'B', 'C', 'D',])
        self.domain_taxa.lock()
        self.range_taxa.lock()
        self.mapping_func = lambda x: self.range_taxa.require_taxon(label=x.label[0].upper())
        self.mapping_dict = {}
        for t in self.domain_taxa:
            self.mapping_dict[t] = self.mapping_func(t)
            t.containing_taxa = self.mapping_dict[t]
        self.expected_forward_label_map = {
                'a1' : 'A',
                'a2' : 'A',
                'a3' : 'A',
                'a4' : 'A',
                'a5' : 'A',
                'a6' : 'A',
                'b1' : 'B',
                'b2' : 'B',
                'b3' : 'B',
                'b4' : 'B',
                'c1' : 'C',
                'c2' : 'C',
                'c3' : 'C',
                'c4' : 'C',
                'd1' : 'D',
                'd2' : 'D',
                'd3' : 'D',}
        self.expected_backward_label_map = {
                'A' : set(['a1', 'a2', 'a3', 'a4', 'a5', 'a6']),
                'B' : set(['b1', 'b2', 'b3', 'b4',]),
                'C' : set(['c1', 'c2', 'c2', 'c3',]),
                'D' : set(['d1', 'd2', 'd3'])
                }

    def verifyMapping(self, tsm):
        for t in self.domain_taxa:
            self.assertEqual(tsm.forward[t].label, self.expected_forward_label_map[t.label])
        for t in self.range_taxa:
            self.assertEqual(set([i.label for i in tsm.reverse[t]]), self.expected_backward_label_map[t.label])

    def testFromFunc(self):
        tsm = dendropy.TaxonSetMapping(mapping_func=self.mapping_func, domain_taxon_set=self.domain_taxa)
        self.verifyMapping(tsm)

    def testFromAttr(self):
        tsm = dendropy.TaxonSetMapping(mapping_attr_name='containing_taxa', domain_taxon_set=self.domain_taxa)
        self.verifyMapping(tsm)

    def testFromDict(self):
        tsm = dendropy.TaxonSetMapping(mapping_dict=self.mapping_dict)
        self.verifyMapping(tsm)

if __name__ == "__main__":
    unittest.main()


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
        self.assertDistinctButEqual(ts, self.taxon_set)

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

class TaxonSetPartitionTest(ExtendedTestCase):

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
            self.membership_lists = [
            [self.taxon_set[0], self.taxon_set[1], self.taxon_set[2], self.taxon_set[3],
             self.taxon_set[13], self.taxon_set[14]],
            [self.taxon_set[4], self.taxon_set[5], self.taxon_set[6], self.taxon_set[7]],
            [self.taxon_set[8], self.taxon_set[9], self.taxon_set[10], self.taxon_set[11]],
            [self.taxon_set[12], self.taxon_set[15], self.taxon_set[16]]
        ]
        self.list_index_to_label_map = ['a', 'b', 'c', 'd']

    def verify_membership_func(self, tsp_mfunc):
        for t in self.taxon_set:
            self.assertEqual(self.membership_func(t), tsp_mfunc(t))

    def verify_membership_dict(self, tsp_md):
        for k, v in tsp_md.items():
            self.assertIn(k, self.membership_dict)
            self.assertEqual(self.membership_dict[k], v)

    def verify_membership_lists(self, tsp_ml):
        self.assertEqual(len(self.membership_lists), len(tsp_ml))
        for i in tsp_ml:
            self.assertIn(i, self.membership_lists)

    def testFromMembershipFunc(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_func=self.membership_func)
        self.verify_membership_func(tsp.membership_func)
        self.verify_membership_dict(tsp.get_membership_dict())
        self.verify_membership_lists(tsp.get_membership_lists())

    def testFromMembershipDict(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_dict=self.membership_dict)
        self.verify_membership_func(tsp.membership_func)
        self.verify_membership_dict(tsp.get_membership_dict())
        self.verify_membership_lists(tsp.get_membership_lists())

    def testFromMembershipLists(self):
        tsp = dendropy.TaxonSetPartition(self.taxon_set, membership_lists=self.membership_lists)
        tsp_mfunc = tsp.membership_func
        for t in self.taxon_set:
            self.assertEqual(self.membership_func(t), self.list_index_to_label_map[tsp_mfunc(t)])
        tsp_md = tsp.get_membership_dict()
        for k, v in tsp_md.items():
            self.assertIn(k, self.membership_dict)
            self.assertEqual(self.membership_dict[k], self.list_index_to_label_map[v])
        self.verify_membership_lists(tsp.get_membership_lists())

if __name__ == "__main__":
    unittest.main()

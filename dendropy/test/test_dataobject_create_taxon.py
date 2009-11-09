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
import dendropy

class TaxaTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.labels = []
        for idx in xrange(10):
            self.labels.append("T%d" % (idx+1))
        self.taxon_set = dendropy.TaxonSet()
        for label in self.labels:
            self.taxon_set.new_taxon(label=label)

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
        self.assertSame(x1, self.taxon_set.get_taxon(label="X1"))
        self.assertSame(self.taxon_set.get_taxon(label="X2"), None)
        self.taxon_set.require_taxon(label="X3")
        self.assertEquals(len(self.taxon_set), 12)

    def testTaxonQuerying(self):
        ts = dendropy.TaxonSet(self.labels)
        self.assertSame(ts.get_taxon(label="Q"), None)
        self.assertSame(ts.get_taxon(label="T1"), ts[0])

if __name__ == "__main__":
    unittest.main()

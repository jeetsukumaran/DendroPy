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
from dendropy.test.support import framework
import dendropy

class TaxaTest(framework.DataObjectVerificationTestCase):

    def setUp(self):
        self.labels = []
        for idx in xrange(10):
            self.labels.append("T%d" % (idx+1))
        self.taxon_set = dendropy.TaxonSet()
        for label in self.labels:
            self.taxon_set.new_taxon(label=label)

    def testLocking(self):
        self.taxon_set.lock()
        self.assertEquals(len(self.taxon_set), 10)
        for idx, t in enumerate(self.taxon_set):
            self.assertEquals(t.label, self.labels[idx])
        self.assertRaises(KeyError, self.taxon_set.new_taxon, label="A1")
        self.assertRaises(KeyError, self.taxon_set.require_taxon, label="A1", oid=None)
        self.taxon_set.unlock()
        x1 = self.taxon_set.new_taxon(label="X1")
        self.assertIsSame(x1, self.taxon_set.get_taxon(label="X1"))
        self.assertIsSame(self.taxon_set.get_taxon(label="X2"), None)
        self.taxon_set.require_taxon(label="X3")
        self.assertEquals(len(self.taxon_set), 12)

    def testCompositionFromStrings(self):
        ts = dendropy.TaxonSet(self.labels)
        self.assertDistinctButEqual(ts, self.taxon_set)

    def x(self):
        self.assertTrue(ts.has_taxa(labels=self.labels))
        self.assertFalse(ts.has_taxa(labels=labels+["Z"]))
        self.assertTrue(ts.has_taxa(taxa=ts))
        self.assertFalse(ts.has_taxa(taxa=tj))
        self.assertFalse(ts.has_taxon(label="Z"))
        self.assertFalse(ts.has_taxon(taxon=Taxon()))
        for label in labels:
            self.assertTrue(ts.has_taxon(label=label))

    def testFromListOfStrings(self):
        pass

if __name__ == "__main__":
    unittest.main()

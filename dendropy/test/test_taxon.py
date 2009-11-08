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
Tests creation, addition, and deletion of taxa.
"""

import unittest
from dendropy.utility import messaging
from dendropy.dataobject import taxon

_LOG = messaging.get_logger(__name__)

class TaxaTest(unittest.TestCase):

    def runTest(self):
        _LOG.info("Testing TaxonSet")
        ti = taxon.TaxonSet()
        labels = []
        for idx in xrange(10):
            labels.append("T%d" % (idx+1))
        for label in labels:
            ti.new_taxon(label=label)
        self.assertEquals(len(ti), 10)
        _LOG.info(ti)
        for idx, t in enumerate(ti):
            self.assertEquals(t.label, labels[idx])
        ti.lock()
        self.assertRaises(Exception, ti.new_taxon, label="A1")
        self.assertRaises(Exception, ti.require_taxon, label="A1", oid=None)
        ti.unlock()
        ti.new_taxon("X1")
        self.assertEquals(ti.get_taxon(label="X2"), None)
        ti.require_taxon(label="X3")
        self.assertEquals(len(ti), 12)

        tj = taxon.TaxonSet()
        tj1 = tj.new_taxon(label="Z")

        self.assertTrue(ti.has_taxa(labels=labels))
        self.assertFalse(ti.has_taxa(labels=labels+["Z"]))
        self.assertTrue(ti.has_taxa(taxa=ti))
        self.assertFalse(ti.has_taxa(taxa=tj))
        self.assertFalse(ti.has_taxon(label="Z"))
        self.assertFalse(ti.has_taxon(taxon=tj1))
        for label in labels:
            self.assertTrue(ti.has_taxon(label=label))

if __name__ == "__main__":
    unittest.main()

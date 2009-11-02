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

    def test_taxonobj(self):
        _LOG.info("Testing TaxonSet")
        ti = taxon.TaxonSet()
        for idx in xrange(10):
            ti.new_taxon("T%d" % (idx+1))
        self.assertEquals(len(ti), 10)
        _LOG.info(ti)
        for idx, t in enumerate(ti):
#             _LOG.info("'%s'=='T%d'?" % (t.label, idx+1))
            self.assertEquals(t.label, "T%d" % (idx+1))
        ti.lock()
        self.assertRaises(Exception, ti.new_taxon, label="A1")
        self.assertRaises(Exception, ti.require_taxon, label="A1", oid=None)
        ti.unlock()
        ti.new_taxon("X1")
        self.assertEquals(ti.get_taxon("X2"), None)
        ti.require_taxon("X3")
        self.assertEquals(len(ti), 12)

if __name__ == "__main__":
    unittest.main()

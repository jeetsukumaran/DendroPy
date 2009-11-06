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
Verifies that data objects generated for use in testing are correct.
"""

import unittest
from dendropy.test.support import datagen
from dendropy.test.support import framework
import dendropy

class TreeInstantiationTest(framework.DataObjectVerificationTestCase):

    def testTreeFromStandard(self):
        tree1 = datagen.get_standard_four_taxon_tree()
        node_oids = [nd.oid for nd in tree1.postorder_node_iter()]
        self.assertEqual(node_oids, ['a', 'b', 'i1', 'c', 'd', 'i2', 'root'])
        tax_labels = [nd.taxon.label for nd in tree1.postorder_node_iter() if nd.taxon is not None]
        self.assertEqual(tax_labels, ['A', 'B', 'C', 'D'])

if __name__ == "__main__":
    unittest.main()

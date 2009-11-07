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
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import framework
import dendropy

class DataForTestingTest(framework.DataObjectVerificationTestCase):

    def testTreeFromStandard(self):
        tree1 = datagen.four_taxon_tree1()
        node_oids = [nd.oid for nd in tree1.postorder_node_iter()]
        self.assertEqual(node_oids, ['a', 'b', 'i1', 'c', 'd', 'i2', 'root'])
        tax_labels = [nd.taxon.label for nd in tree1.postorder_node_iter() if nd.taxon is not None]
        self.assertEqual(tax_labels, ['A', 'B', 'C', 'D'])

    def testReferenceTreeList(self):
        tlist1 = datagen.reference_tree_list()
        ref_trees_newick = [n.strip() for n in datagen.reference_tree_list_newick_string().split(";")]
        ref_node_labels = datagen.reference_tree_list_postorder_node_labels()
        ref_node_rels = datagen.reference_tree_list_node_relationships()
        for ti, t1 in enumerate(tlist1):
            t1.assign_node_labels_from_taxon_or_oid()

            t1_newick = t1.as_newick_str(include_internal_labels=True)
            self.assertEqual(t1_newick, ref_trees_newick[ti])

            node_labels1 = [nd.label for nd in t1.postorder_node_iter()]
            self.assertEqual(node_labels1, ref_node_labels[ti])

            nodes1 = [nd for nd in t1.postorder_node_iter()]
            for ndi, nd1 in enumerate(nodes1):
                ndrel = ref_node_rels[ti][nd1.label]
                ndrel.test_node(self, nd1)

if __name__ == "__main__":
    unittest.main()

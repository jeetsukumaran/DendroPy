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
from dendropy.test.support import datatest
import dendropy

class DataForTestingTest(datatest.DataObjectVerificationTestCase):

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

    def testReferenceDnaArray(self):
        char_array = datagen.reference_dna_array()
        taxon_set = char_array.taxon_set
        dna_dict = datagen.reference_dna_dict()
        self.assertEqual(len(char_array), 13)
        self.assertSame(char_array.default_state_alphabet, dendropy.DNA_STATE_ALPHABET)
        for tax_label, tax_seq_symbols in dna_dict.items():
            taxon = taxon_set.require_taxon(label=tax_label)
            self.assertContained(taxon, char_array)
            seq_vec = char_array[taxon]
            self.assertEqual(len(seq_vec), len(tax_seq_symbols))
            for si, s1 in enumerate(tax_seq_symbols):
                state = seq_vec[si].value
                self.assertEqual(state.symbol, s1)
            self.assertEqual(seq_vec.symbols_as_string(), tax_seq_symbols)

    def testReferenceSingleTaxonSetDataSet(self):
        dataset = datagen.reference_single_taxonset_dataset()
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.taxon_sets[0]), 13)
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertSame(dataset.tree_lists[0].taxon_set, dataset.taxon_sets[0])
        self.assertEqual(len(dataset.char_arrays), 2)
        self.assertSame(dataset.char_arrays[0].taxon_set, dataset.taxon_sets[0])
        self.assertSame(dataset.char_arrays[1].taxon_set, dataset.taxon_sets[0])

if __name__ == "__main__":
    unittest.main()

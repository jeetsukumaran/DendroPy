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
            t1_newick = t1.as_newick_string(include_internal_labels=True, preserve_spaces=True)
            self.assertEqual(t1_newick, ref_trees_newick[ti])
            node_labels1 = [nd.label for nd in t1.postorder_node_iter()]
            self.assertEqual(node_labels1, ref_node_labels[ti])
            nodes1 = [nd for nd in t1.postorder_node_iter()]
            for ndi, nd1 in enumerate(nodes1):
                ndrel = ref_node_rels[ti][nd1.label]
                ndrel.test_node(self, nd1)

    def testReferenceDnaMatrix(self):
        char_matrix = datagen.reference_dna_matrix()
        taxon_set = char_matrix.taxon_set
        dna_dict = datagen.reference_dna_dict()
        self.assertEqual(len(char_matrix), 33)
        self.assertIs(char_matrix.default_state_alphabet, dendropy.DNA_STATE_ALPHABET)
        for tax_label, tax_seq_symbols in dna_dict.items():
            taxon = taxon_set.require_taxon(label=tax_label)
            self.assertIn(taxon, char_matrix)
            seq_vec = char_matrix[taxon]
            self.assertEqual(len(seq_vec), len(tax_seq_symbols))
            for si, s1 in enumerate(tax_seq_symbols):
                state = seq_vec[si].value
                self.assertEqual(state.symbol, s1)
            self.assertEqual(seq_vec.symbols_as_string(), tax_seq_symbols)

    def testReferenceSingleTaxonSetDataSet(self):
        dataset = datagen.reference_single_taxonset_dataset()
        self.assertEqual(len(dataset.taxon_sets), 1)
        self.assertEqual(len(dataset.taxon_sets[0]), 33)
        self.assertEqual(len(dataset.tree_lists), 1)
        self.assertIs(dataset.tree_lists[0].taxon_set, dataset.taxon_sets[0])
        self.assertEqual(len(dataset.char_matrices), 2)
        self.assertIs(dataset.char_matrices[0].taxon_set, dataset.taxon_sets[0])
        self.assertIs(dataset.char_matrices[1].taxon_set, dataset.taxon_sets[0])

if __name__ == "__main__":
    unittest.main()

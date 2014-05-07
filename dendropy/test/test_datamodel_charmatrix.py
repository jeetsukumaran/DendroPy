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
Tests character sequence map.
"""

import unittest
import dendropy
from dendropy.datamodel import charmatrixmodel
from dendropy.test.support import dendropytest

class TaxonCharacterMatrixBasicCrud(dendropytest.ExtendedTestCase):

    def setUp(self):
        self.taxon_namespace = dendropy.TaxonNamespace()
        for i in range(3):
            label = "T{}".format(i)
            t = self.taxon_namespace.require_taxon(label=label)

    def test_setitem(self):
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=self.taxon_namespace)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_seq_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(self.taxon_namespace)
        for idx, taxon in enumerate(self.taxon_namespace):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_seq_map), len(self.taxon_namespace))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_seq_map))
        for idx, taxon in enumerate(self.taxon_namespace):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

if __name__ == "__main__":
    unittest.main()

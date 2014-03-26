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
Tests for NEWICK taxon handling.
"""

import sys
import os
import unittest
import dendropy
from dendropy.dataio import nexusprocessing

class TaxonSymbolMappingTest(unittest.TestCase):

    def test_standard_lookup_and_create(self):
        labels = ["t{}".format(i) for i in range(1, 10)]
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        for idx, label in enumerate(labels):
            self.assertEqual(len(tns), idx)
            t1 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertEqual(t1.label, label)
            t2 = tsm.require_taxon_for_symbol(label)
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t2)
            self.assertEqual(t2.label, label)
            t3 = tsm.require_taxon_for_symbol(str(idx+1))
            self.assertEqual(len(tns), idx+1)
            self.assertIs(t1, t3)
            self.assertEqual(t3.label, label)

    def test_taxon_namespace_locking(self):
        tns = dendropy.TaxonNamespace()
        tsm = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=tns)
        self.assertFalse(tns.is_mutable)
        del tsm
        self.assertTrue(tns.is_mutable)


if __name__ == "__main__":
    unittest.main()

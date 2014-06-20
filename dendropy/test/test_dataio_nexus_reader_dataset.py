# !/usr/bin/env python

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
NEXUS data read/write parse/format tests.
"""

from dendropy.test.support import pathmap
from dendropy.test.support import dendropytest
from dendropy.utility import messaging
import unittest
import dendropy
_LOG = messaging.get_logger(__name__)

class DataSetNexusTaxonManagement(dendropytest.ExtendedTestCase):

    def testMultiTaxonSet(self):
        d = dendropy.DataSet()
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[1]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 3)
        self.assertEqual(len(d.taxon_sets[2]), 33)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 4)
        self.assertEqual(len(d.taxon_sets[3]), 114)

    def testBoundTaxonSetDefault(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonSet()
        d.attach_taxon_set(t)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertIs(d.taxon_sets[0], d.attached_taxon_set)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 147)

    def testBindAndUnbind(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonSet()
        d.attach_taxon_set(t)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertIs(d.taxon_sets[0], d.attached_taxon_set)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        _LOG.info(d.taxon_sets[0].description(2))
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.detach_taxon_set()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        self.assertEqual(len(d.taxon_sets[1]), 114)

    def testBindToSpecifiedTaxonSet(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonSet()
        d.attach_taxon_set(t)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertIs(d.taxon_sets[0], d.attached_taxon_set)
        self.assertIs(d.attached_taxon_set, t)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        d.detach_taxon_set()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[0]), 33)
        self.assertEqual(len(d.taxon_sets[1]), 114)

if __name__ == "__main__":
    unittest.main()

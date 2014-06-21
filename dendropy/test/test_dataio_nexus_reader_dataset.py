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

    def testMultiTaxonNamespace(self):
        d = dendropy.DataSet()
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 2)
        self.assertEqual(len(d.taxon_namespaces[1]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_namespaces), 3)
        self.assertEqual(len(d.taxon_namespaces[2]), 33)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_namespaces), 4)
        self.assertEqual(len(d.taxon_namespaces[3]), 114)

    def testBoundTaxonNamespaceDefault(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonNamespace()
        d.attach_taxon_namespace(t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 147)

    def testBindAndUnbind(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonNamespace()
        d.attach_taxon_namespace(t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        _LOG.info(d.taxon_namespaces[0].description(2))
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.detach_taxon_namespace()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_namespaces), 2)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        self.assertEqual(len(d.taxon_namespaces[1]), 114)

    def testBindToSpecifiedTaxonNamespace(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonNamespace()
        d.attach_taxon_namespace(t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        self.assertIs(d.attached_taxon_namespace, t)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.detach_taxon_namespace()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_namespaces), 2)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        self.assertEqual(len(d.taxon_namespaces[1]), 114)

if __name__ == "__main__":
    unittest.main()

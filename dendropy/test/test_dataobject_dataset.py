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
Tests creation, reading, update, deletion of DataSet objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
from dendropy.test.support import pathmap
from dendropy.test.support import datatest
from dendropy.test.support import datagen
import dendropy

class DataSetCreateTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.dataset = datagen.reference_single_taxonset_dataset()

    def testReadNexus(self):
        ds_str = self.dataset.as_string(schema="nexus")
        ds2 = dendropy.DataSet(stream=StringIO(ds_str), schema="nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromFileFactory(self):
        ds_str = self.dataset.as_string(schema="nexus")
        ds2 = dendropy.DataSet.get_from_stream(StringIO(ds_str), "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromPathFactory(self):
        ds2 = dendropy.DataSet.get_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromStringFactory(self):
        ds_str = self.dataset.as_string(schema="nexus")
        ds2 = dendropy.DataSet.get_from_string(ds_str, "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromCopy(self):
        ds2 = dendropy.DataSet(self.dataset)
        self.assertDistinctButEqual(self.dataset, ds2)

    def testSimpleCopyDnaMatrix(self):
        char_matrix = datagen.reference_dna_matrix()
        ds1 = dendropy.DataSet(char_matrix)
        self.assertEqual(len(ds1.char_matrices), 1)
        self.assertIs(ds1.char_matrices[0], char_matrix)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

    def testSimpleCopyStandardMatrix(self):
        char_matrix = datagen.reference_standard_matrix()
        ds1 = dendropy.DataSet(char_matrix)
        self.assertEqual(len(ds1.char_matrices), 1)
        self.assertIs(ds1.char_matrices[0], char_matrix)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

class DataSetTaxonManagement(datatest.DataObjectVerificationTestCase):

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
        d = dendropy.DataSet(attach_taxon_set=True)
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
        d = dendropy.DataSet(attach_taxon_set=True)
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

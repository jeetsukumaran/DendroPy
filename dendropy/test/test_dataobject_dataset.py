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
Tests creation, reading, update, deletion of DataSet objects.
"""

import unittest
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.test.support import pathmap
from dendropy.test.support import datatest
from dendropy.test.support import datagen
import dendropy

class DataSetCreateTest(datatest.DataObjectVerificationTestCase):

    def setUp(self):
        self.dataset = datagen.reference_single_taxonset_dataset()

    def testReadNexus(self):
        ds_str = self.dataset.as_string(format="nexus")
        ds2 = dendropy.DataSet(stream=StringIO(ds_str), format="nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromFileFactory(self):
        ds_str = self.dataset.as_string(format="nexus")
        ds2 = dendropy.DataSet.get_from_stream(StringIO(ds_str), "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromPathFactory(self):
        ds2 = dendropy.DataSet.get_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromStringFactory(self):
        ds_str = self.dataset.as_string(format="nexus")
        ds2 = dendropy.DataSet.get_from_string(ds_str, "nexus")
        self.assertDistinctButEqual(self.dataset, ds2)

    def testFromCopy(self):
        ds2 = dendropy.DataSet(self.dataset)
        self.assertDistinctButEqual(self.dataset, ds2)

    def testSimpleCopyDnaMatrix(self):
        char_matrix = datagen.reference_dna_matrix()
        ds1 = dendropy.DataSet(char_matrix)
        self.assertEqual(len(ds1.char_matrices), 1)
        self.assertSame(ds1.char_matrices[0], char_matrix)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

    def testSimpleCopyStandardMatrix(self):
        char_matrix = datagen.reference_standard_matrix()
        ds1 = dendropy.DataSet(char_matrix)
        self.assertEqual(len(ds1.char_matrices), 1)
        self.assertSame(ds1.char_matrices[0], char_matrix)
        ds2 = dendropy.DataSet(ds1)
        self.assertDistinctButEqual(ds1, ds2)

class DataSetTaxonManagement(datatest.DataObjectVerificationTestCase):

    def testMultiTaxonSet(self):
        d = dendropy.DataSet(multi_taxon_set=True)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[1]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 3)
        self.assertEqual(len(d.taxon_sets[2]), 29)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 4)
        self.assertEqual(len(d.taxon_sets[3]), 114)

    def testBoundTaxonSetDefault(self):
        d = dendropy.DataSet(attached_taxon_set=True)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertSame(d.taxon_sets[0], d.attached_taxon_set)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 143)

    def testBindAndUnbind(self):
        d = dendropy.DataSet(attached_taxon_set=True)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertSame(d.taxon_sets[0], d.attached_taxon_set)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.detach_taxon_set()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        self.assertEqual(len(d.taxon_sets[1]), 114)

    def testBindToSpecifiedTaxonSet(self):
        d = dendropy.DataSet(multi_taxon_set=True)
        t = dendropy.TaxonSet()
        d.attach_taxon_set(t)
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertSame(d.taxon_sets[0], d.attached_taxon_set)
        self.assertSame(d.attached_taxon_set, t)
        d.read_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.read_from_path(pathmap.tree_source_path('pythonidae.reference-trees.newick'), "newick")
        self.assertEqual(len(d.taxon_sets), 1)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        d.detach_taxon_set()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), "proteinfasta")
        self.assertEqual(len(d.taxon_sets), 2)
        self.assertEqual(len(d.taxon_sets[0]), 29)
        self.assertEqual(len(d.taxon_sets[1]), 114)

if __name__ == "__main__":
    unittest.main()

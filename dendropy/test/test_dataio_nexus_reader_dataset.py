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
from dendropy.test.support import standard_file_test_chars
from dendropy.utility import messaging
import unittest
import dendropy
_LOG = messaging.get_logger(__name__)

class DataSetNexusCharsGetFromTestCase(dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.check_taxon_annotations = False
        cls.check_matrix_annotations = False
        cls.check_sequence_annotations = False
        cls.check_column_annotations = False
        cls.check_cell_annotations = False
        standard_file_test_chars.DnaTestChecker.build()
        standard_file_test_chars.RnaTestChecker.build()
        standard_file_test_chars.ProteinTestChecker.build()
        standard_file_test_chars.Standard01234TestChecker.build()

    def test_single(self):
        srcs = (
                ("standard-test-chars-dna.multi.nexus", standard_file_test_chars.DnaTestChecker),
                ("standard-test-chars-rna.multi.nexus", standard_file_test_chars.RnaTestChecker),
                ("standard-test-chars-protein.multi.nexus", standard_file_test_chars.ProteinTestChecker),
                ("standard-test-chars-generic.interleaved.nexus", standard_file_test_chars.Standard01234TestChecker),
                )
        for src_filename, src_matrix_checker_type in srcs:
            src_path = pathmap.char_source_path(src_filename)
            ds = dendropy.DataSet.get_from_path(src_path,
                    "nexus")
            self.assertEqual(len(ds.char_matrices), 1)
            self.assertEqual(len(ds.taxon_namespaces), 1)
            self.assertIs(ds.char_matrices[0].taxon_namespace,
                    ds.taxon_namespaces[0])
            self.assertEqual(type(ds.char_matrices[0]), src_matrix_checker_type.matrix_type)
            if src_matrix_checker_type.matrix_type is dendropy.StandardCharacterMatrix:
                src_matrix_checker_type.create_class_data_label_sequence_map_based_on_state_alphabet(src_matrix_checker_type,
                        ds.char_matrices[0].default_state_alphabet)
            standard_file_test_chars.general_char_matrix_checker(
                    self,
                    ds.char_matrices[0],
                    src_matrix_checker_type,
                    check_taxon_annotations=self.check_taxon_annotations,
                    check_matrix_annotations=self.check_matrix_annotations,
                    check_sequence_annotations=self.check_sequence_annotations,
                    check_column_annotations=self.check_column_annotations,
                    check_cell_annotations=self.check_cell_annotations,)

class DataSetNexusTaxonManagementTestCase(dendropytest.ExtendedTestCase):

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

    def testAttachTaxonNamespaceOnGet(self):
        t = dendropy.TaxonNamespace()
        d = dendropy.DataSet.get_from_path(pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'),
                "nexus",
                taxon_namespace=t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIsNot(d.attached_taxon_namespace, None)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        self.assertIs(d.attached_taxon_namespace, t)
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

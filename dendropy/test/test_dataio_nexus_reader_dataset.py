# !/usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
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
from dendropy.test.support import curated_test_tree
from dendropy.test.support import standard_file_test_datasets
from dendropy.utility import messaging
import unittest
import dendropy
_LOG = messaging.get_logger(__name__)

class DataSetNexusSingleCharsTestCase(dendropytest.ExtendedTestCase):

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
        standard_file_test_chars.ContinuousTestChecker.build()
        cls.srcs = (
                ("standard-test-chars-dna.multi.nexus", standard_file_test_chars.DnaTestChecker),
                ("standard-test-chars-rna.multi.nexus", standard_file_test_chars.RnaTestChecker),
                ("standard-test-chars-protein.multi.nexus", standard_file_test_chars.ProteinTestChecker),
                ("standard-test-chars-generic.interleaved.nexus", standard_file_test_chars.Standard01234TestChecker),
                ("standard-test-chars-continuous.mesquite.nexus", standard_file_test_chars.ContinuousTestChecker),
                )

    def verify_char_matrix(self, char_matrix, src_matrix_checker_type):
        self.assertEqual(type(char_matrix), src_matrix_checker_type.matrix_type)
        if src_matrix_checker_type.matrix_type is dendropy.StandardCharacterMatrix:
            src_matrix_checker_type.create_class_fixtures_label_sequence_map_based_on_state_alphabet(src_matrix_checker_type,
                    char_matrix.default_state_alphabet)
        standard_file_test_chars.general_char_matrix_checker(
                self,
                char_matrix,
                src_matrix_checker_type,
                check_taxon_annotations=self.check_taxon_annotations,
                check_matrix_annotations=self.check_matrix_annotations,
                check_sequence_annotations=self.check_sequence_annotations,
                check_column_annotations=self.check_column_annotations,
                check_cell_annotations=self.check_cell_annotations,)

    def test_get_single(self):
        for src_filename, src_matrix_checker_type in self.__class__.srcs:
            src_path = pathmap.char_source_path(src_filename)
            ds = dendropy.DataSet.get_from_path(src_path, "nexus")
            self.assertEqual(len(ds.char_matrices), 1)
            self.assertEqual(len(ds.taxon_namespaces), 1)
            self.assertIs(ds.char_matrices[0].taxon_namespace,
                    ds.taxon_namespaces[0])
            self.verify_char_matrix(ds.char_matrices[0], src_matrix_checker_type)

    def test_read_single(self):
        for src_filename, src_matrix_checker_type in self.__class__.srcs:
            src_path = pathmap.char_source_path(src_filename)
            ds = dendropy.DataSet()
            result = ds.read(path=src_path, schema="nexus")
            self.assertEqual(result, (1,0,1))
            self.assertEqual(len(ds.char_matrices), 1)
            self.assertEqual(len(ds.taxon_namespaces), 1)
            self.assertIs(ds.char_matrices[0].taxon_namespace,
                    ds.taxon_namespaces[0])
            self.verify_char_matrix(ds.char_matrices[0], src_matrix_checker_type)

    def test_read_successive_unattached_taxon_namespace(self):
        ds = dendropy.DataSet()
        for src_idx, (src_filename, src_matrix_checker_type) in enumerate(self.__class__.srcs):
            src_path = pathmap.char_source_path(src_filename)
            result = ds.read(path=src_path, schema="nexus")
            self.assertEqual(result, (1,0,1))
            self.assertEqual(len(ds.char_matrices), src_idx+1)
            self.assertEqual(len(ds.taxon_namespaces), src_idx+1)
            self.assertIs(ds.char_matrices[src_idx].taxon_namespace,
                    ds.taxon_namespaces[src_idx])
            self.verify_char_matrix(ds.char_matrices[src_idx], src_matrix_checker_type)

    def test_read_successive_attached_taxon_namespace(self):
        ds = dendropy.DataSet()
        tns = dendropy.TaxonNamespace()
        ds.attach_taxon_namespace(tns)
        for src_idx, (src_filename, src_matrix_checker_type) in enumerate(self.__class__.srcs):
            src_path = pathmap.char_source_path(src_filename)
            result = ds.read(path=src_path, schema="nexus")
            self.assertEqual(result, (0,0,1))
            self.assertEqual(len(ds.char_matrices), src_idx+1)
            self.assertEqual(len(ds.taxon_namespaces), 1)
            self.assertIs(ds.char_matrices[src_idx].taxon_namespace,
                    ds.taxon_namespaces[0])
            self.assertIs(ds.taxon_namespaces[0], tns)
            self.verify_char_matrix(ds.char_matrices[src_idx], src_matrix_checker_type)

class DataSetNexusMultipleCharBlocksTestCase(dendropytest.ExtendedTestCase):

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

    def verify_dataset(self, ds):
        self.assertEqual(len(ds.taxon_namespaces), 1)
        tns = ds.taxon_namespaces[0]
        checkers = (
                standard_file_test_chars.RnaTestChecker,
                standard_file_test_chars.ProteinTestChecker,
                standard_file_test_chars.Standard01234TestChecker,
                standard_file_test_chars.DnaTestChecker,
                )
        self.assertEqual(len(ds.char_matrices), len(checkers))
        for idx, (char_matrix, checker) in enumerate(zip(ds.char_matrices, checkers)):
            self.assertIs(char_matrix.taxon_namespace, tns)
            if checker.matrix_type is dendropy.StandardCharacterMatrix:
                checker.create_class_fixtures_label_sequence_map_based_on_state_alphabet(checker,
                        char_matrix.default_state_alphabet)
            standard_file_test_chars.general_char_matrix_checker(
                    self,
                    char_matrix,
                    checker,
                    check_taxon_annotations=self.check_taxon_annotations,
                    check_matrix_annotations=self.check_matrix_annotations,
                    check_sequence_annotations=self.check_sequence_annotations,
                    check_column_annotations=self.check_column_annotations,
                    check_cell_annotations=self.check_cell_annotations,)

    def test_basic_get(self):
        src_filename = "standard-test-chars-multiple-char-blocks.1.basic.nexus"
        src_path = pathmap.char_source_path(src_filename)
        ds = dendropy.DataSet.get_from_path(src_path, "nexus")
        self.verify_dataset(ds)

    def test_basic_read(self):
        src_filename = "standard-test-chars-multiple-char-blocks.1.basic.nexus"
        src_path = pathmap.char_source_path(src_filename)
        ds = dendropy.DataSet()
        result = ds.read(path=src_path, schema="nexus")
        self.assertEqual(result, (1, 0, 4))
        self.verify_dataset(ds)

class DataSetNexusTreesTestCase(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def test_multiple_trees(self):
        src_filename = "multitreeblocks.nex"
        src_path = pathmap.tree_source_path(src_filename)
        ds = dendropy.DataSet.get_from_path(src_path, "nexus")
        self.assertEqual(len(ds.taxon_namespaces), 1)
        self.assertEqual(len(ds.tree_lists), 3)

    def test_with_translate(self):
        srcs = (
            ("curated-with-translate-block-and-internal-taxa.nex", False),
            ("curated-with-translate-block-and-untranslated-internal-taxa.nex", True),
            )
        for src_filename, suppress_internal_taxa in srcs:
            src_path = pathmap.tree_source_path(src_filename)
            ds = dendropy.DataSet.get_from_path(
                    src_path,
                    "nexus",
                    suppress_internal_node_taxa=suppress_internal_taxa)
            self.assertEqual(len(ds.tree_lists), 1)
            tree_list = ds.tree_lists[0]
            tree_labels = ("1", "2", "3")
            self.assertEqual(len(tree_list), len(tree_labels))
            for tree_idx, (tree, label) in enumerate(zip(tree_list, tree_labels)):
                self.assertEqual(tree.label, label)
                self.verify_curated_tree(
                        tree=tree,
                        suppress_internal_node_taxa=suppress_internal_taxa,
                        suppress_leaf_node_taxa=False,
                        suppress_edge_lengths=False,
                        node_taxon_label_map=None)

class DataSetNexusReaderMixedTestCase(
        standard_file_test_datasets.StandardSingleTaxonNamespaceDataSet,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_datasets.StandardSingleTaxonNamespaceDataSet.build(cls)
        cls.check_taxon_annotations = False
        cls.check_matrix_annotations = False
        cls.check_sequence_annotations = False
        cls.check_column_annotations = False
        cls.check_cell_annotations = False

    def test_basic_get(self):
        src_filename = "standard-test-mixed.1.basic.nexus"
        src_path = pathmap.mixed_source_path(src_filename)
        ds = dendropy.DataSet.get_from_path(
                src_path,
                "nexus",
                suppress_internal_node_taxa=False, # so internal labels get translated
                )
        self.verify_dataset(ds)

    def test_basic_read(self):
        src_filename = "standard-test-mixed.1.basic.nexus"
        src_path = pathmap.mixed_source_path(src_filename)
        ds = dendropy.DataSet()
        result = ds.read(
                path=src_path,
                schema="nexus",
                suppress_internal_node_taxa=False, # so internal labels get translated
                )
        self.assertEqual(result, (1, 7, 4))
        self.verify_dataset(ds)

class DataSetNexusTaxonManagementTestCase(dendropytest.ExtendedTestCase):

    def testMultiTaxonNamespace(self):
        d = dendropy.DataSet()
        d.read(path=pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.mle.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 2)
        self.assertEqual(len(d.taxon_namespaces[1]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.reference-trees.newick'), schema="newick")
        self.assertEqual(len(d.taxon_namespaces), 3)
        self.assertEqual(len(d.taxon_namespaces[2]), 33)
        d.read(path=pathmap.char_source_path('caenophidia_mos.chars.fasta'), schema="fasta", data_type="protein")
        self.assertEqual(len(d.taxon_namespaces), 4)
        self.assertEqual(len(d.taxon_namespaces[3]), 114)

    def testBoundTaxonNamespaceDefault(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonNamespace()
        d.attach_taxon_namespace(t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        d.read(path=pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.mle.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.reference-trees.newick'), schema="newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.char_source_path('caenophidia_mos.chars.fasta'), schema="fasta", data_type="protein")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 147)

    def testBindAndUnbind(self):
        d = dendropy.DataSet()
        t = dendropy.TaxonNamespace()
        d.attach_taxon_namespace(t)
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertIs(d.taxon_namespaces[0], d.attached_taxon_namespace)
        self.assertIs(d.attached_taxon_namespace, t)
        d.read(path=pathmap.mixed_source_path('reference_single_taxonset_dataset.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.mle.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.reference-trees.newick'), schema="newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.detach_taxon_namespace()
        d.read(path=pathmap.char_source_path('caenophidia_mos.chars.fasta'), schema="fasta", data_type="protein")
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
        d.read(path=pathmap.tree_source_path('pythonidae.mle.nex'), schema="nexus")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.read(path=pathmap.tree_source_path('pythonidae.reference-trees.newick'), schema="newick")
        self.assertEqual(len(d.taxon_namespaces), 1)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        d.detach_taxon_namespace()
        d.read_from_path(pathmap.char_source_path('caenophidia_mos.chars.fasta'), schema="fasta", data_type="protein")
        self.assertEqual(len(d.taxon_namespaces), 2)
        self.assertEqual(len(d.taxon_namespaces[0]), 33)
        self.assertEqual(len(d.taxon_namespaces[1]), 114)

class NexusDataSetMultiTreesBlockTestCase(dendropytest.ExtendedTestCase):

    def test_multiple_trees1(self):
        src_filename = "multitreeblocks.nex"
        src_path = pathmap.tree_source_path(src_filename)
        ds = dendropy.DataSet.get_from_path(src_path, "nexus")
        self.assertEqual(len(ds.tree_lists), 3)
        for tt in ds.tree_lists:
            self.assertEqual(len(tt), 3)

    def test_multiple_trees2(self):
        src_filename = "multitreeblocks2.nex"
        src_path = pathmap.tree_source_path(src_filename)
        ds = dendropy.DataSet.get_from_path(src_path, "nexus")
        self.assertEqual(len(ds.taxon_namespaces), 1)
        self.assertEqual(len(ds.tree_lists), 2)
        for tt in ds.tree_lists:
            self.assertEqual(len(tt), 2)
            self.assertIs(tt.taxon_namespace, ds.taxon_namespaces[0])

class DataSetNexusReaderMesquiteMultipleTaxonNamespacesTest(
        standard_file_test_datasets.MultipleTaxonNamespaceDataSet,
        dendropytest.ExtendedTestCase):

    def test_attached_taxon_namespace(self):
        tns = dendropy.TaxonNamespace()
        ds = dendropy.DataSet.get_from_path(
                pathmap.mixed_source_path('multitaxa_mesquite.nex'),
                "nexus",
                taxon_namespace=tns)
        self.verify_attached_taxon_namespace(ds, tns)

    def test_default(self):
        ds = dendropy.DataSet.get_from_path(
                pathmap.mixed_source_path('multitaxa_mesquite.nex'),
                "nexus")
        self.verify_unrestricted(ds)

if __name__ == "__main__":
    unittest.main()

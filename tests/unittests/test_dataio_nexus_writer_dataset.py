#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Tests for NEXUS dataset writing.
"""

import unittest
import dendropy
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap
from support import dendropytest
from support import compare_and_validate
from support import standard_file_test_datasets

class DataSetNexusWriterMixedTestCase(
        standard_file_test_datasets.StandardSingleTaxonNamespaceDataSet,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_datasets.StandardSingleTaxonNamespaceDataSet.build(cls)
        cls.check_taxon_annotations = False
        cls.check_matrix_annotations = False
        cls.check_sequence_annotations = False
        cls.check_column_annotations = False
        cls.check_cell_annotations = False

    def get_source_dataset(self, **kwargs):
        src_ds = dendropy.DataSet.get_from_path(
                pathmap.mixed_source_path("standard-test-mixed.1.basic.nexus"),
                "nexus",
                **kwargs
                )
        return src_ds

    def test_basic(self):
        # `suppress_internal_node_taxa=False` so internal labels get translated
        d0 = self.get_source_dataset(suppress_internal_node_taxa=False)
        s = self.write_out_validate_equal_and_return(
                d0, "nexus", {})
        ds = dendropy.DataSet.get_from_string(
                s, "nexus",
                suppress_internal_node_taxa=False, # so internal labels get translated
                )
        self.verify_dataset(ds)

class DataSetNexusWriterMesquiteMultipleTaxonNamespacesTest(
        standard_file_test_datasets.MultipleTaxonNamespaceDataSet,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    def test_attached_taxon_namespace(self):
        d0 = dendropy.DataSet.get_from_path(
                pathmap.mixed_source_path('multitaxa_mesquite.nex'),
                "nexus")
        for tns in d0.taxon_namespaces:
            d0.attach_taxon_namespace(tns)
            s = self.write_out_validate_equal_and_return(
                    d0, "nexus", {})
            ds = dendropy.DataSet.get_from_string(s, "nexus",)
            self.verify_attached_taxon_namespace_written(ds, tns)

    def test_default(self):
        d0 = dendropy.DataSet.get_from_path(
                pathmap.mixed_source_path('multitaxa_mesquite.nex'),
                "nexus")
        s = self.write_out_validate_equal_and_return(
                d0, "nexus", {})
        ds = dendropy.DataSet.get_from_string(
                s, "nexus",
                )
        self.verify_unrestricted(ds)

if __name__ == "__main__":
    unittest.main()

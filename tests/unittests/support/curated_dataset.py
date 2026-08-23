#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
DataSet test data generation and verification.
"""

import dendropy
from . import standard_file_test_chars

class CuratedDataSetGenerator(object):

    def get_dataset(self):
        ds = dendropy.DataSet()
        tns1 = dendropy.TaxonNamespace(label="tns1")

    def get_taxon_namespace(self):
        pass

    def get_tree_list(self, taxon_namespace):
        pass

    def get_standard_char_matrix(self, taxon_namespace):
        pass

    def get_dna_char_matrix(self, taxon_namespace):
        pass

    def get_rna_char_matrix(self, taxon_namespace):
        pass

    def get_protein_char_matrix(self, taxon_namespace):
        pass



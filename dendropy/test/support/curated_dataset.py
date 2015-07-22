#! /usr/bin/env python

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
DataSet test data generation and verification.
"""

import dendropy
from dendropy.test.support import standard_file_test_chars

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



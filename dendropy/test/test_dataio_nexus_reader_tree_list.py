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
Tests for general NEXUS tree list reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NexusStandardTreeListReaderTestCase(
        standard_file_test_trees.StandardTreeListReaderTestCase,
        dendropytest.ExtendedTestCase):

    @classmethod
    def build(cls):
        standard_file_test_trees.StandardTreeListReaderTestCase.build(schema="nexus",
                taxa_on_tree_equal_taxa_in_taxon_namespace=False)

    @classmethod
    def setUpClass(cls):
        cls.build()

    ## NOTE: tests are in standard_file_test_trees.StandardTreeListReaderTestCase !! ##

if __name__ == "__main__":
    unittest.main()

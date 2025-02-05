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
Tests for general NEWICK tree iteration reading.
"""

import sys
import unittest
import dendropy
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest
from support import standard_file_test_trees

class NewickTreeIteratorReaderDefaultTestCase(
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls)

if __name__ == "__main__":
    unittest.main()

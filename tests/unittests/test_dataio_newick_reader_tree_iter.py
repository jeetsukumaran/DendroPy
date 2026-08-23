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

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
Tests for D3 writing.
"""

import collections
import unittest
import re
import dendropy
from dendropy.utility.textprocessing import StringIO
from dendropy.test.support import pathmap

class D3WriterSingleTreeTest(unittest.TestCase):

    def test_basic(self):
        tree = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path('pythonidae.mle.nex'), schema="nexus")
        s = StringIO()
        d3_tree_str = tree.write(
                file=s,
                schema="d3")

class D3WriterTreeListTest(unittest.TestCase):

    def test_basic(self):
        trees = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path("pythonidae.reference-trees.newick"), "newick")
        s = StringIO()
        d3_trees_str = trees.write(
                file=s,
                schema="d3")


if __name__ == "__main__":
    unittest.main()


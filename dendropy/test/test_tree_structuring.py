#! /usr/bin/env python

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
Tests native tree structuring routines.
"""

import unittest
from dendropy.test.support import pathmap
from dendropy.utility import messaging
import dendropy

_LOG = messaging.get_logger(__name__)

class TestTreeRestructures(unittest.TestCase):

    def testLadderizeLeft(self):
        """
        Dummy test!
        """
        tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        tree.ladderize()

    def testLadderizeRight(self):
        """
        Dummy test!
        """
        tree = dendropy.Tree.get_from_path(pathmap.tree_source_path('pythonidae.mle.nex'), "nexus")
        tree.ladderize(right=True)

if __name__ == "__main__":
    unittest.main()


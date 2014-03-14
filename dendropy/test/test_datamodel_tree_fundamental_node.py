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
Tests basic Node operations.
"""

import unittest
import dendropy

class TestNode(unittest.TestCase):

    def test_basic_construction(self):
        taxon = dendropy.Taxon("z")
        nd = dendropy.Node(taxon=taxon, label="x", edge_length=1)
        self.assertIs(nd.taxon, taxon)
        self.assertEqual(nd.label, "x")
        edge = nd.edge
        self.assertEqual(edge.length, 1)
        self.assertIs(edge.head_node, nd)

if __name__ == "__main__":
    unittest.main()

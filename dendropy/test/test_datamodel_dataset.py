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
Tests basic DataSet curation
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest

class DataSetCreation(dendropytest.ExtendedTestCase):

    def test_basic_create_and_add(self):
        ds = dendropy.DataSet(label="d1")
        self.assertEqual(ds.label, label)

if __name__ == "__main__":
    unittest.main()


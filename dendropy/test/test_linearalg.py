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
Tests linear algebra and matrix computation.
"""

import unittest
from dendropy.test.support import extendedtest
from dendropy.mathlib import linearalg
from dendropy.utility import messaging

_LOG = messaging.get_logger(__name__)

class TestMatrix(unittest.TestCase):

    def setUp(self):
        self.g1 = [
            [2, 2],
            [2, 5],
            [6, 5],
            [7, 3],
            [4, 7],
            [6, 4],
            [5, 3],
            [4, 6],
            [2, 5],
            [1, 3],
            ]
        self.g2 = [
            [6, 5],
            [7, 4],
            [8, 7],
            [5, 6],
            [5, 4],
            ]
        self.h1 = list(linearalg.new_matrix(self.g1).tr())
        self.h2 = list(linearalg.new_matrix(self.g2).tr())
        self.cov1 = linearalg.new_matrix([ [3.89, 0.13],[0.13, 2.21] ])
        self.cov2 = linearalg.new_matrix([ [1.36, 0.56],[0.56, 1.36] ])
        self.expected_cov = [self.cov1, self.cov2]

    def testCovarianceByColumns(self):
        for i, x in enumerate([self.g1, self.g2]):
            x = linearalg.new_matrix(x)
            s = x.covariance_by_cols(population_variance=True)
            assert self.expected_cov[i] == s

    def testCovarianceByRows(self):
        for i, x in enumerate([self.h1, self.h2]):
            x = linearalg.new_matrix(x)
            s = x.covariance_by_rows(population_variance=True)
            assert self.expected_cov[i] == s

if __name__ == "__main__":
    unittest.main()

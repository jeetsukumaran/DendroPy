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
from dendropy.mathlib import distance
from dendropy.mathlib import linearalg
from dendropy.utility import messaging

_LOG = messaging.get_logger(__name__)

class TestMahalanobis(unittest.TestCase):

    def testMahalanobisDistance2d(self):
        g1 = [
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
        g2 = [
            [6, 5],
            [7, 4],
            [8, 7],
            [5, 6],
            [5, 4],
            ]
        d = distance.mahalanobis(g1, g2, population_variance=True)
        self.assertAlmostEqual(d, 1.410417840)

    def testMahalanobisDistance2d_MATLAB(self):
        # MATLAB code:
        #    g1 = [ 2 2; 2 5; 6 5; 7 3; 4 7; 6 4; 5 3; 4 6; 2 5; 1 3; ];
        #    g2 = [ 6 5; 7 4; 8 7; 5 6; 5 4; ];
        #    mahal(g1, g2)
        g1 = [ [2, 2], [2, 5], [6, 5], [7, 3], [4, 7], [6, 4], [5, 3], [4, 6], [2, 5], [1, 3], ]
        g2 = [ [6, 5], [7, 4], [8, 7], [5, 6], [5, 4], ]
        e =  [ 11.9083, 12.0333, 0.0333, 4.9083, 8.0333, 0.9083, 2.9083, 4.9083, 12.0333, 15.9083 ]
        s = linearalg.new_matrix(g2).covariance_by_cols(population_variance=False)
        for i, g in enumerate(g1):
            self.assertAlmostEqual(e[i], distance.squared_mahalanobis([g], g2, cov=s), 4)

if __name__ == "__main__":
    unittest.main()

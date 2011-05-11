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
Tests statistical routines.
"""

import unittest
from dendropy.test.support import extendedtest
from dendropy.utility import statistics
from dendropy.utility import messaging

_LOG = messaging.get_logger(__name__)

class HpdCalculationTests(extendedtest.ExtendedTestCase):

    def testEmpiricalUnimodalHpd1(self):
        v = [100.93356, 100.66576, 99.44097, 100.60761, 103.65723, 101.15563,
                99.09657, 100.39654, 98.77339, 101.13712, 99.33979, 99.99060,
                100.39395, 101.68240, 100.99664, 99.17798, 100.83020, 98.90373,
                100.30441, 99.49553, 100.52652, 99.76291, 99.95605, 99.63605,
                99.21535, 100.51619, 100.55036, 101.21747, 101.04181, 97.76084,
                100.19069, 99.46182, 100.47579, 99.56889, 100.23977, 101.22907,
                97.85931, 100.86051, 99.56121, 100.44109, 100.02328, 98.62446,
                100.11008, 100.12700, 99.27087, 100.72895, 99.06796, 99.38019,
                99.79908, 100.82761, 101.26901, 99.88911, 98.09761, 99.16706,
                98.98752, 100.10088, 100.58883, 99.42982, 101.90322, 101.22817,
                101.36052, 97.70629, 100.15950, 99.39458, 100.19414, 103.43317,
                100.32429, 98.90429, 101.28049, 99.82948, 100.96041, 99.46024,
                98.22509, 101.63878, 100.66998, 101.82238, 99.49847, 100.41055,
                98.71792, 99.66001, 98.53177, 99.11997, 100.14802, 98.96423,
                101.93145, 100.09478, 100.85930, 99.82181, 101.50284, 99.93301,
                99.57168, 98.19978, 100.90708, 99.25086, 101.74170, 99.86034,
                99.85785, 99.89154, 99.62313, 99.41994,]
        #==> 101.82238 97.70629
        #==> 101.90322 97.76084
        #==> 101.93145 97.85931
        #==> 103.43317 98.09761
        #==> 103.65723 98.19978
        #n= 100
        #nn= 5
        #xx= [4.11609, 4.142380000000003, 4.0721400000000045, 5.335560000000001, 5.457449999999994]
        #m= 4.07214
        #nnn= 100
        #(97.85931, 101.93145)
        c1, c2 = statistics.calc_unimodal_emp_hpd(v)
        self.assertAlmostEqual(c1, 97.85931)
        self.assertAlmostEqual(c2, 101.93145)


if __name__ == "__main__":
    unittest.main()


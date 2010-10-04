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
Test of data collection classes.
"""

import unittest
from dendropy.utility import containers
from dendropy.test.support.extendedtest import ExtendedTestCase

class TestNormalizedBitmaskDict(ExtendedTestCase):

    def runTest(self):
        """Testing NormalizedBitmaskDict"""
        mask = 0xFF # 1111 1111
        splits = [
            ((0x03, '0000 0011'), (0x03, '0000 0011')),
            ((0x34, '0011 0100'), (0xCB, '1100 1011')),
            ((0x44, '0100 0100'), (0xBB, '1011 1011')),
            ((0x12, '0001 0010'), (0xED, '1110 1101')),
            ((0x75, '0111 0101'), (0x75, '0111 0101')),
            ]
        d = containers.NormalizedBitmaskDict(mask=mask)
        for s in splits:
            d[s[0][0]] = s[0][1]

        for s in splits:
            self.assertIn(s[0][0], d)
            self.assertIn(s[1][0], d)
            self.assertEqual(d[s[0][0]], d[s[1][0]])

        for k, v in d.items():
            pass

        del d[splits[0][0][0]]
        del d[splits[1][1][0]]
        self.assertNotIn(splits[0][0][0], d)
        self.assertNotIn(splits[0][1][0], d)
        self.assertNotIn(splits[1][0][0], d)
        self.assertNotIn(splits[1][1][0], d)

class TestOrderedSet(unittest.TestCase):

    class DummyObject(object):

        def __init__(self, label):
            self.label = label

        def __str__(self):
            return "<%s: '%s'>" % (id(self), self.label)

    def runTest(self):
        """Testing OrderedSet [DUMMY TEST -- TODO!]"""
        d = containers.OrderedSet()
        for i in xrange(10):
            d.add(TestOrderedSet.DummyObject(i))
        ### TODO! ###

if __name__ == "__main__":
    unittest.main()

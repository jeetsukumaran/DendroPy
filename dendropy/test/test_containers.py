#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

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
            self.assertContained(s[0][0], d)
            self.assertContained(s[1][0], d)
            self.assertEqual(d[s[0][0]], d[s[1][0]])

        for k, v in d.items():
            pass

        del d[splits[0][0][0]]
        del d[splits[1][1][0]]
        self.assertNotContained(splits[0][0][0], d)
        self.assertNotContained(splits[0][1][0], d)
        self.assertNotContained(splits[1][0][0], d)
        self.assertNotContained(splits[1][1][0], d)

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

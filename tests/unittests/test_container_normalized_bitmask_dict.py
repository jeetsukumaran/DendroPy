#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
NormalizedBitmaskDict tests.
"""

import unittest
from dendropy.utility import container

class TestNormalizedBitmaskDict(unittest.TestCase):

    def runTest(self):
        """Testing NormalizedBitmaskDict"""
        fill_bitmask = 0xFF # 1111 1111
        splits = [
            ((0x03, '0000 0011'), (0x03, '0000 0011')),
            ((0x34, '0011 0100'), (0xCB, '1100 1011')),
            ((0x44, '0100 0100'), (0xBB, '1011 1011')),
            ((0x12, '0001 0010'), (0xED, '1110 1101')),
            ((0x75, '0111 0101'), (0x75, '0111 0101')),
            ]
        d = container.NormalizedBitmaskDict(fill_bitmask=fill_bitmask)
        for s in splits:
            d[s[0][0]] = s[0][1]

        for s in splits:
            self.assertIn(s[0][0], d)
            self.assertIn(s[1][0], d)
            self.assertEqual(d[s[0][0]], d[s[1][0]])

        # Regression test: NormalizedBitmaskDict.__setitem__() used to call
        # dict.__setitem__() directly, bypassing OrderedDict's internal
        # bookkeeping. This left len()/__repr__() correct but silently
        # broke iteration (.items()/.keys()/.values()/__iter__()), which
        # would all report zero entries despite the dict being non-empty.
        seen_keys = set(k for k, v in d.items())
        self.assertEqual(len(seen_keys), len(splits))
        self.assertEqual(len(d), len(splits))

        del d[splits[0][0][0]]
        del d[splits[1][1][0]]
        self.assertNotIn(splits[0][0][0], d)
        self.assertNotIn(splits[0][1][0], d)
        self.assertNotIn(splits[1][0][0], d)
        self.assertNotIn(splits[1][1][0], d)

        # Regression test: __delitem__() had the same dict-vs-OrderedDict
        # bypass problem; confirm iteration still reflects the deletions.
        seen_keys_after_delete = set(k for k, v in d.items())
        self.assertEqual(len(seen_keys_after_delete), len(splits) - 2)
        self.assertEqual(len(d), len(splits) - 2)

if __name__ == "__main__":
    unittest.main()

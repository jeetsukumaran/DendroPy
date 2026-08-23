#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Tests for the ``dendropy.application.sumtrees`` command-line application.
"""

import unittest

from dendropy.application import sumtrees


class AvailableCpuCountTestCase(unittest.TestCase):

    def test_returns_positive_int(self):
        count = sumtrees._available_cpu_count()
        self.assertIsInstance(count, int)
        self.assertGreaterEqual(count, 1)


if __name__ == "__main__":
    unittest.main()

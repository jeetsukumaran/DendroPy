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
Discrete character tests.
"""

import random
import unittest
from dendropy.model import discrete
from dendropy.simulate import treesim

class DiscreteCharacterEvolverTest(unittest.TestCase):

    def test_rng_param(self):
        """Test determinism from rng param.

        Regression test for #170.
        """
        tree = treesim.birth_death_tree(
            birth_rate=1.0,
            death_rate=0.5,
            num_extant_tips=4,
            rng=random.Random(100)
        )
        assert (
            discrete.hky85_chars(10, tree, rng=random.Random(1)).as_string(
                schema="nexml",
            )
            == discrete.hky85_chars(10, tree, rng=random.Random(1)).as_string(
                schema="nexml",
            )
        )
        assert (
            discrete.hky85_chars(10, tree, rng=random.Random(1)).as_string(
                schema="nexml",
            )
            != discrete.hky85_chars(10, tree, rng=random.Random(2)).as_string(
                schema="nexml",
            )
        )

if __name__ == "__main__":
    unittest.main()

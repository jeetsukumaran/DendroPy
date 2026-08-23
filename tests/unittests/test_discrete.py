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
Discrete character tests.
"""

import random
import unittest
import dendropy
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

    def test_simulate_discrete_chars_honors_root_states(self):
        """
        Regression test: simulate_discrete_chars() accepts a 'root_states'
        parameter but was passing 'root_states=None' on to the underlying
        evolver, silently discarding the caller's requested root sequence.
        """
        tree = dendropy.Tree.get(data="(A:1,B:1):0;", schema="newick")
        seq_model = discrete.Jc69()
        alphabet_states = list(seq_model.state_alphabet)
        seq_len = 8
        root_states = [alphabet_states[i % 4] for i in range(seq_len)]
        discrete.simulate_discrete_chars(
                seq_len=seq_len,
                tree_model=tree,
                seq_model=seq_model,
                root_states=root_states,
                retain_sequences_on_tree=True,
                rng=random.Random(1),
                )
        root_node = tree.seed_node
        root_sequence = root_node.sequences[-1]
        self.assertEqual(list(root_sequence), root_states)

if __name__ == "__main__":
    unittest.main()

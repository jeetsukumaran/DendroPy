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
Tests of birth-death model likelihood calculatio and fitting.
"""

import unittest
import itertools as it
import json
import os
import sys
import dendropy
from dendropy.model import birthdeath
sys.path.insert(0, os.path.dirname(__file__))
from support.mockrandom import MockRandom
from support import pathmap

class PureBirthModelEstimationTestCase(unittest.TestCase):

    def testBasicEst(self):
        # list of tuples, (birth-rate, log-likelihood)
        expected_results = (
                # birth rate               # log-likelihood
            (0.02879745490817826186758, -59.41355682054444287132355),
            (0.03074708092192806122012, -57.38280732060526645454956),
            (0.02539588437187430269848, -63.31025321526630023072357),
            (0.02261951969802362960582, -66.89924384677527768872096),
            (0.02804607815688910446572, -60.23314120509648716961237),
            (0.02748663302756114423797, -60.85775993426526042640035),
            (0.02816256618562208019485, -60.10465085978295007862471),
            (0.03592126646048716259729, -52.56123967307649991198559),
            (0.02905144990609926855529, -59.14133401672411594063306),
            (0.02703739196351075124714, -61.36860953277779628933786),
            (0.01981322730236481297061, -71.00561162515919022553135),
        )
        trees = dendropy.TreeList.get_from_path(
                pathmap.tree_source_path("pythonidae.reference-trees.newick"), "newick")
        self.assertEqual(len(trees), len(expected_results))
        for tree, expected_result in zip(trees, expected_results):
            obs_result1 = birthdeath.fit_pure_birth_model(tree=tree, ultrametricity_precision=1e-5)
            obs_result2 = birthdeath.fit_pure_birth_model(internal_node_ages=tree.internal_node_ages(ultrametricity_precision=1e-5))
            for obs_result in (obs_result1, obs_result2):
                self.assertAlmostEqual(obs_result["birth_rate"], expected_result[0], 5)
                self.assertAlmostEqual(obs_result["log_likelihood"], expected_result[1], 5)

class BirthDeathTreeTest(unittest.TestCase):
    def testGSABD(self):
        """test that the birth-death process produces the correct number of tips with GSA."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 15), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.2, num_extant_tips=num_leaves, gsa_ntax=3*num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testYule(self):
        """test that the pure-birth process produces the correct number of tips."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testGSA(self):
        """test that the pure-birth process produces the correct number of tips with GSA."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=num_leaves, gsa_ntax=4*num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testBDTree(self):
        """PureCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.2, num_extant_tips=num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testBDTreeTime(self):
        """Time-based stop condition -- tree generation without checking [TODO: checks]"""
        _RNG = MockRandom()
        for tree_factory in [dendropy.Tree, lambda: None]:
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.2, max_time=4, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())


class FastBirthDeathTreeTest(unittest.TestCase):
    def testGSABD(self):
        """test that the birth-death process produces the correct number of tips with GSA."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 15), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.2, num_extant_tips=num_leaves, gsa_ntax=3*num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testYule(self):
        """test that the pure-birth process produces the correct number of tips."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.fast_birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testGSA(self):
        """test that the pure-birth process produces the correct number of tips with GSA."""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.fast_birth_death_tree(birth_rate=1.0, death_rate=0.0, num_extant_tips=num_leaves, gsa_ntax=4*num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testBDTree(self):
        """PureCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        _RNG = MockRandom()
        for num_leaves, tree_factory in it.product(
            range(2, 20), [dendropy.Tree, lambda: None]
        ):
            t = birthdeath.fast_birth_death_tree(birth_rate=1.0, death_rate=0.2, num_extant_tips=num_leaves, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEqual(num_leaves, len(t.leaf_nodes()))

    def testBDTreeTime(self):
        """Time-based stop condition -- tree generation without checking [TODO: checks]"""
        _RNG = MockRandom()
        for tree_factory in [dendropy.Tree, lambda: None]:
            t = birthdeath.fast_birth_death_tree(birth_rate=1.0, death_rate=0.2, max_time=4, tree=tree_factory(), rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())

class BirthDeathLikelihoodTestCases(unittest.TestCase):

    def test_likelihood_calc(self):
        src = pathmap.other_source_stream("birth-death-test-data1.json")
        ref_data = json.load(src)
        for test_group in ref_data:
            tree = dendropy.Tree.get(
                    data=test_group["tree"] + ";",
                    schema="newick",
                    rooting="force-rooted")
            estimation_profiles = test_group["estimates"]
            for estimation_profile in estimation_profiles:
                expected = estimation_profile["log_likelihood"]
                observed1 = birthdeath.birth_death_likelihood(
                        tree=tree,
                        birth_rate=estimation_profile["estimation_birth_rate"],
                        death_rate=estimation_profile["estimation_death_rate"],
                        sampling_probability=estimation_profile["estimation_sampling_probability"],
                        sampling_strategy=estimation_profile["estimation_sampling_strategy"],
                        is_mrca_included=estimation_profile["estimation_includes_mrca"],
                        condition_on=estimation_profile["estimation_conditioned_on"],
                        ultrametricity_precision=1e-5,
                        )
                self.assertAlmostEqual(observed1, expected, 6)
                observed2 = birthdeath.birth_death_likelihood(
                        internal_node_ages=sorted(tree.internal_node_ages(ultrametricity_precision=1e-5), reverse=True),
                        birth_rate=estimation_profile["estimation_birth_rate"],
                        death_rate=estimation_profile["estimation_death_rate"],
                        sampling_probability=estimation_profile["estimation_sampling_probability"],
                        sampling_strategy=estimation_profile["estimation_sampling_strategy"],
                        is_mrca_included=estimation_profile["estimation_includes_mrca"],
                        condition_on=estimation_profile["estimation_conditioned_on"],
                        )
                self.assertAlmostEqual(observed2, expected, 6)

if __name__ == "__main__":
    unittest.main()


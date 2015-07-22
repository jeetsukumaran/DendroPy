#! /usr/bin/env python

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
Tests state alphabet definition and management.
"""

import sys
import random
import itertools
import unittest
import collections
import dendropy
from dendropy.utility import container
from dendropy.test.support import dendropytest

class StateAlphabetTester(object):

    def validate_state_identities(self,
            state_container,
            state_iter,
            symbol_iter,
            expected_symbols,
            expected_denomination,
            member_state_map,
            additional_synonyms_map,
            case_sensitive=True):
        self.assertEqual(len(state_container), len(expected_symbols))
        states = list(state_iter())
        self.assertEqual(len(states), len(expected_symbols))
        canonical_symbols = list(symbol_iter(include_synonyms=False))
        self.assertEqual(canonical_symbols, expected_symbols)
        for state, symbol in zip(states, expected_symbols):
            all_synonyms = []
            self.assertEqual(state.symbol, symbol)
            self.assertEqual(state.state_denomination, expected_denomination)
            if member_state_map:
                expected_member_state_symbols = tuple(member_state_map[symbol])
                self.assertEqual(state.fundamental_symbols, expected_member_state_symbols)
                fundamental_states = self.sa.get_fundamental_states_for_symbols(state.symbol)
                fss = [fs.symbol for fs in fundamental_states]
                self.assertEqual(tuple(fss), expected_member_state_symbols)
            else:
                self.assertEqual(state.fundamental_states, (state,))
                self.assertEqual(state.fundamental_symbols, (state.symbol,))
            if case_sensitive:
                self.assertNotIn(symbol.upper(), state.symbol_synonyms)
                self.assertNotIn(symbol.lower(), state.symbol_synonyms)
            else:
                if symbol.upper() != symbol:
                    self.assertIn(symbol.upper(), state.symbol_synonyms)
                    all_synonyms.append(symbol.upper())
                if symbol.lower() != symbol:
                    self.assertIn(symbol.lower(), state.symbol_synonyms)
                    all_synonyms.append(symbol.lower())
            if additional_synonyms_map:
                x = additional_synonyms_map.get(state.symbol, None)
                if x is not None:
                    all_synonyms.extend(x)
            self.assertEqual(set(all_synonyms), set(state.symbol_synonyms),
                    state)
            self.assertEqual(len(all_synonyms), len(state.symbol_synonyms), state.symbol)

            expected_fundamental_state_symbols = member_state_map.get(state.symbol, None)
            if expected_fundamental_state_symbols is None:
                expected_fundamental_states = set()
            else:
                expected_fundamental_states = self.sa.get_states_for_symbols(expected_fundamental_state_symbols)
                check_ss = [x.symbol for x in expected_fundamental_states]
                self.assertEqual(set(check_ss), set(expected_fundamental_state_symbols))

    def test_fundamental_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._fundamental_states,
                state_iter=self.sa.fundamental_state_iter,
                symbol_iter=self.sa.fundamental_symbol_iter,
                expected_symbols=self.expected_fundamental_state_symbols,
                expected_denomination=self.sa.FUNDAMENTAL_STATE,
                member_state_map={},
                additional_synonyms_map=self.additional_synonyms_map,
                case_sensitive=False)

    def test_ambiguous_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._ambiguous_states,
                state_iter=self.sa.ambiguous_state_iter,
                symbol_iter=self.sa.ambiguous_symbol_iter,
                expected_symbols=self.expected_ambiguous_state_symbols,
                expected_denomination=self.sa.AMBIGUOUS_STATE,
                member_state_map=self.ambiguous_symbol_mappings,
                additional_synonyms_map=self.additional_synonyms_map,
                case_sensitive=False)

    def test_polymorphic_state_definitions(self):
        self.validate_state_identities(
                state_container=self.sa._polymorphic_states,
                state_iter=self.sa.polymorphic_state_iter,
                symbol_iter=self.sa.polymorphic_symbol_iter,
                expected_symbols=self.expected_polymorphic_state_symbols,
                expected_denomination=self.sa.POLYMORPHIC_STATE,
                member_state_map=self.polymorphic_symbol_mappings,
                additional_synonyms_map=self.additional_synonyms_map,
                case_sensitive=False)

    def test_state_iter(self):
        states = list(self.sa.state_iter())
        self.assertEqual(len(states), self.num_total_states)
        self.assertEqual(len(self.sa), len(states))
        expected_state_symbol_iter = itertools.chain(
                self.expected_fundamental_state_symbols,
                self.ambiguous_symbol_mappings,
                self.polymorphic_symbol_mappings
                )
        for state, symbol in zip(states, expected_state_symbol_iter):
            self.assertEqual(state.symbol, symbol)

    def test_symbol_iter(self):
        # assumes that the state iterators -- fundamental_state_iter,
        # ambiguous_state_iter, etc. -- all work as advertised
        iter_groups = (
                (self.sa.fundamental_symbol_iter, self.sa.fundamental_state_iter),
                (self.sa.ambiguous_symbol_iter, self.sa.ambiguous_state_iter),
                (self.sa.polymorphic_symbol_iter, self.sa.polymorphic_state_iter),
                (self.sa.multistate_symbol_iter, self.sa.multistate_state_iter),
                )
        for symbol_iter, state_iter in iter_groups:
            states = list(state_iter())
            for include_synonyms in (False, True):
                expected_symbols = []
                for state in states:
                    if state.symbol:
                        expected_symbols.append(state.symbol)
                    if include_synonyms:
                        for ss in state.symbol_synonyms:
                            expected_symbols.append(ss)
                obs_symbols = list(symbol_iter(include_synonyms=include_synonyms))
                self.assertEqual(expected_symbols, obs_symbols)

    def test_symbol_state_pair_iter(self):
        states = list(self.sa.state_iter())
        for include_synonyms in (False, True):
            expected_pairs = []
            for state in states:
                if state.symbol:
                    expected_pairs.append((state.symbol, state,))
                if include_synonyms:
                    for ss in state.symbol_synonyms:
                        expected_pairs.append((ss, state,))
            obs_pairs = list(self.sa.symbol_state_pair_iter(include_synonyms=include_synonyms))
            self.assertEqual(expected_pairs, obs_pairs)

    def test_state_denomination(self):
        for state in self.sa.fundamental_state_iter():
            self.assertEqual(state.state_denomination, self.sa.FUNDAMENTAL_STATE)
        for state in self.sa.ambiguous_state_iter():
            self.assertEqual(state.state_denomination, self.sa.AMBIGUOUS_STATE)
        for state in self.sa.polymorphic_state_iter():
            self.assertEqual(state.state_denomination, self.sa.POLYMORPHIC_STATE)

    def test_compiled_lookup_immutability(self):
        self.sa.compile_lookup_mappings()
        for m in (
                self.sa.canonical_symbol_state_map,
                self.sa.full_symbol_state_map,
                self.sa._fundamental_states_to_ambiguous_state_map,
                self.sa._fundamental_states_to_polymorphic_state_map
                ):
            if m:
                k = list(m.keys())[0]
            else:
                k = 1
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                m[k] = 1
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                del m[k]
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                m.pop(k)
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                m.clear()
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                m.update({})
            with self.assertRaises(container.FrozenOrderedDict.ImmutableTypeError):
                m.fromkeys([1,2,3])
        # check if re-compilation is possible
        self.sa.compile_lookup_mappings()

    def test_canonical_symbol_state_map(self):
        m = self.sa.canonical_symbol_state_map
        states = list(self.sa.state_iter())
        exp_symbols = [s.symbol for s in states if s.symbol]
        obs_symbols = list(m)
        self.assertEqual(obs_symbols, exp_symbols)
        self.assertEqual(len(m), len(states))
        for obs_symbol, exp_state in zip(m, states):
            self.assertEqual(obs_symbol, exp_state.symbol)
            self.assertIs(m[obs_symbol], exp_state)

    def test_full_symbol_state_map(self):
        m = self.sa.full_symbol_state_map
        states = list(self.sa.state_iter())
        exp_symbols = []
        exp_symbol_state_pairs = []
        for state in states:
            if state.symbol:
                exp_symbols.append(state.symbol)
                exp_symbol_state_pairs.append((state.symbol, state))
                for s in state.symbol_synonyms:
                    exp_symbols.append(s)
                    exp_symbol_state_pairs.append((s, state))
        obs_symbols = list(m)
        if self.sa.no_data_state is not None:
            exp_symbols.insert(0, None)
            exp_symbol_state_pairs.insert(0, (None, self.sa.no_data_state))
        self.assertEqual(obs_symbols, exp_symbols)
        self.assertEqual(len(m), len(exp_symbols))
        self.assertEqual(len(m), len(exp_symbol_state_pairs))
        for obs_symbol, exp_symbol, sspair in zip(m, exp_symbols, exp_symbol_state_pairs):
            self.assertEqual(obs_symbol, exp_symbol)
            self.assertEqual(obs_symbol, sspair[0])
            self.assertIs(m[obs_symbol], sspair[1])

    def test_no_data_state(self):
        if self.sa.no_data_state is not None:

            # some setup
            expected_fundamental_states = list(self.sa.fundamental_state_iter())
            expected_fundamental_symbols = [s.symbol for s in expected_fundamental_states]
            test_symbols = [None] + expected_fundamental_symbols
            expected_states = [self.sa.no_data_state] + expected_fundamental_states

            # check definitions
            self.assertIn(self.sa.no_data_state, self.sa._ambiguous_states)
            self.assertEqual(self.sa.no_data_state.symbol, self.sa.no_data_symbol)
            self.assertEqual(self.sa.no_data_state._member_states, tuple(expected_fundamental_states))

            # check look-up map
            full_map = self.sa.full_symbol_state_map
            self.assertIn(None, full_map)
            self.assertIs(full_map[None], self.sa.no_data_state)

            # __getitem__
            self.assertIs(self.sa[None], self.sa.no_data_state)

            # get_states_for_symbols
            s = self.sa.get_states_for_symbols(test_symbols)
            self.assertEqual(s, expected_states)
            self.assertIs(s[0], self.sa.no_data_state)

            # get_fundamental_states_for_symbols
            self.assertEqual(self.sa.get_fundamental_states_for_symbols([None]), expected_fundamental_states)

            # get_canonical_symbol_for_symbol
            self.assertEqual(self.sa.get_canonical_symbol_for_symbol(None), self.sa.no_data_symbol)

            # match_ambiguous_state
            self.assertIs(self.sa.match_ambiguous_state(expected_fundamental_symbols), self.sa.no_data_state)

        else:
            full_map = self.sa.full_symbol_state_map
            self.assertNotIn(None, full_map)

    def test_getitem(self):
        alphabet = self.sa
        for state in self.sa.state_iter():
            self.assertIs(alphabet[state.symbol], state)
            for ss in state.symbol_synonyms:
                self.assertIs(alphabet[ss], state)
            if state._index is not None:
                self.assertIs(alphabet[state._index], state)

#     def test_get_states_for_symbol(self):
#         states = list(self.sa.state_iter())
#         for rep in range(3):
#             n = random.randint(5, 100)
#             selected_states = [self.rng.choice(states) for _ in range(n)]
#             selected_symbols = [s.symbol for s in selected_states]
#             obs_states = self.sa.get_states_for_symbols(selected_symbols)
#             self.assertEqual(obs_states, selected_states)

    def test_get_states_for_symbols(self):
        all_symbols = list(self.sa.full_symbol_state_map.keys())
        for rep in range(3):
            n = random.randint(5, 100)
            selected_symbols = [self.rng.choice(all_symbols) for _ in range(n)]
            selected_states = [self.sa[s] for s in selected_symbols]
            obs_states = self.sa.get_states_for_symbols(selected_symbols)
            self.assertEqual(obs_states, selected_states, "random seed: {}".format(self.random_seed))

    def test_states_property(self):
        check = list(self.sa.state_iter())
        self.assertEqual(len(check), len(self.sa.states))
        for s1, s2 in zip(check, self.sa.states):
            self.assertIs(s1, s2)

    def test_canonical_symbols_property(self):
        check = list(self.sa.canonical_symbol_state_map.keys())
        self.assertEqual(len(check), len(self.sa.states))
        for s1, s2 in zip(check, self.sa.symbols):
            self.assertEqual(s1, s2)

    def test_get_canonical_symbol_for_symbol(self):
        states = list(self.sa.state_iter())
        expected = {}
        no_data_state = None
        for state in states:
            if state.symbol:
                expected[state.symbol] = state.symbol
            if state is self.sa.no_data_state:
                no_data_state = state
            for ss in state.symbol_synonyms:
                expected[ss] = state.symbol
        for symbol in self.sa.full_symbol_state_map:
            if symbol is None:
                self.assertIsNot(self.sa.no_data_state, None)
                self.assertIsNot(self.sa.no_data_symbol, None)
                self.assertIs(self.sa.no_data_state, no_data_state)
            else:
                self.assertEqual(self.sa.get_canonical_symbol_for_symbol(symbol), expected[symbol])

    def test_get_fundamental_states_for_symbols(self):
        all_symbols = list(self.sa.full_symbol_state_map.keys())
        for rep in range(3):
            n = random.randint(5, 100)
            selected_symbols = [self.rng.choice(all_symbols) for _ in range(n)]
            selected_states = []
            for symbol in selected_symbols:
                state = self.sa[symbol]
                if state.state_denomination == self.sa.FUNDAMENTAL_STATE:
                    selected_states.append(state)
                else:
                    if state.state_denomination == self.sa.AMBIGUOUS_STATE:
                        mapping_src = self.ambiguous_symbol_mappings
                    elif state.state_denomination == StateAlphabet.POLYMORPHIC_STATE:
                        mapping_src = self.polymorphic_symbol_mappings
                    else:
                        raise Exception("Unrecognized denomination: {}".format(state.state_denomination))
                    member_states = []
                    canonical_symbol = self.sa.get_canonical_symbol_for_symbol(symbol)
                    for member_symbol in mapping_src[canonical_symbol]:
                        member_states.append(self.sa[member_symbol])
                    selected_states.extend(member_states)
            obs_states = self.sa.get_fundamental_states_for_symbols(selected_symbols)
            if obs_states != selected_states:
                print("\nSelected Symbols: {}\n Selected States: {}\nObserved Symbols: {}\nrandom seed: {}".format(
                        "".join(selected_symbols),
                        "".join([s.symbol for s in selected_states]),
                        "".join([s.symbol for s in obs_states]),
                        self.random_seed))
            self.assertEqual(obs_states, selected_states)

    def test_match_state(self):
        multistate_states = [list(self.sa.ambiguous_state_iter()), list(self.sa.polymorphic_state_iter())]
        match_fns = [self.sa.match_ambiguous_state, self.sa.match_polymorphic_state]
        for multistate_states, match_fn in zip(multistate_states, match_fns):
            for multistate in multistate_states:
                member_states = list(multistate.member_states)
                potential_symbols = []
                for member_state in member_states:
                    member_symbols = [member_state.symbol]
                    for ss in member_state.symbol_synonyms:
                        member_symbols.append(ss)
                    potential_symbols.append(member_symbols)
                for rep in range(5):
                    selected_symbols = [self.rng.choice(x) for x in potential_symbols]
                    self.rng.shuffle(selected_symbols)
                    if self.rng.uniform(0, 1) < 0.5:
                        selected_symbols = "".join(selected_symbols)
                    matched_state = match_fn(selected_symbols)
                    self.assertIs(matched_state, multistate, "random seed: {}".format(self.random_seed))

    def test_on_the_fly_creation_of_multistate(self):
        multistate_states = [list(self.sa.ambiguous_state_iter()), list(self.sa.polymorphic_state_iter())]
        match_fns = [self.sa.match_ambiguous_state, self.sa.match_polymorphic_state]
        add_fns = [self.sa.new_ambiguous_state, self.sa.new_polymorphic_state]
        state_collections = [self.sa._ambiguous_states, self.sa._polymorphic_states]
        symbol_pool = list(self.sa.fundamental_symbol_iter())
        for multistate_states, match_fn, add_fn, state_collection in zip(multistate_states, match_fns, add_fns, state_collections):
            pre_existing_symbol_combinations = []
            new_symbol_combinations = []
            nreps = 0
            while len(new_symbol_combinations) < 3 and nreps < 5:
                nreps += 1
                max_sample_size = min(5, len(self.sa))
                n = self.rng.randint(2, max_sample_size)
                selected_symbols = self.rng.sample(symbol_pool, n)
                try:
                    matched_state = match_fn(selected_symbols)
                except KeyError:
                    new_symbol_combinations.append(selected_symbols)
                    new_state = add_fn(symbol=None, member_state_symbols=selected_symbols)
                    # self.sa.compile_lookup_mappings()
                    try:
                        m2 = match_fn(selected_symbols)
                    except KeyError:
                        raise
                    else:
                        self.assertIs(m2, new_state)
                        self.assertIn(new_state, state_collection, "random seed: {}".format(self.random_seed))
                    finally:
                        state_collection.remove(new_state)
                        self.sa.compile_lookup_mappings()
                else:
                    pre_existing_symbol_combinations.append(selected_symbols)

class DnaStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.random_seed = random.randint(0, sys.maxsize)
        self.rng = random.Random(self.random_seed)

        self.expected_fundamental_state_symbols = ["A", "C", "G", "T", "-"]
        self.ambiguous_symbol_mappings = collections.OrderedDict()
        self.ambiguous_symbol_mappings["?"] = "ACGT-"
        self.ambiguous_symbol_mappings["N"] = "ACGT"
        self.ambiguous_symbol_mappings["R"] = "AG"
        self.ambiguous_symbol_mappings["Y"] = "CT"
        self.ambiguous_symbol_mappings["M"] = "AC"
        self.ambiguous_symbol_mappings["W"] = "AT"
        self.ambiguous_symbol_mappings["S"] = "CG"
        self.ambiguous_symbol_mappings["K"] = "GT"
        self.ambiguous_symbol_mappings["V"] = "ACG"
        self.ambiguous_symbol_mappings["H"] = "ACT"
        self.ambiguous_symbol_mappings["D"] = "AGT"
        self.ambiguous_symbol_mappings["B"] = "CGT"

        self.polymorphic_symbol_mappings = collections.OrderedDict()

        # note reverse polarity here: from referenced to referencing
        self.additional_synonyms_map = collections.OrderedDict()
        self.additional_synonyms_map["N"] = "X"

        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())

        self.sa = dendropy.DNA_STATE_ALPHABET
        self.num_total_states = len(self.expected_fundamental_state_symbols) + len(self.ambiguous_symbol_mappings) + len(self.polymorphic_symbol_mappings)

class RnaStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.random_seed = random.randint(0, sys.maxsize)
        self.rng = random.Random(self.random_seed)

        self.expected_fundamental_state_symbols = ["A", "C", "G", "U", "-"]
        self.ambiguous_symbol_mappings = collections.OrderedDict()
        self.ambiguous_symbol_mappings["?"] = "ACGU-"
        self.ambiguous_symbol_mappings["N"] = "ACGU"
        self.ambiguous_symbol_mappings["R"] = "AG"
        self.ambiguous_symbol_mappings["Y"] = "CU"
        self.ambiguous_symbol_mappings["M"] = "AC"
        self.ambiguous_symbol_mappings["W"] = "AU"
        self.ambiguous_symbol_mappings["S"] = "CG"
        self.ambiguous_symbol_mappings["K"] = "GU"
        self.ambiguous_symbol_mappings["V"] = "ACG"
        self.ambiguous_symbol_mappings["H"] = "ACU"
        self.ambiguous_symbol_mappings["D"] = "AGU"
        self.ambiguous_symbol_mappings["B"] = "CGU"

        self.polymorphic_symbol_mappings = collections.OrderedDict()

        # note reverse polarity here: from referenced to referencing
        self.additional_synonyms_map = collections.OrderedDict()
        self.additional_synonyms_map["N"] = "X"

        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())

        self.sa = dendropy.RNA_STATE_ALPHABET
        self.num_total_states = len(self.expected_fundamental_state_symbols) + len(self.ambiguous_symbol_mappings) + len(self.polymorphic_symbol_mappings)

class NucleotideStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.random_seed = random.randint(0, sys.maxsize)
        self.rng = random.Random(self.random_seed)

        self.expected_fundamental_state_symbols = ["A", "C", "G", "T", "U", "-"]
        self.ambiguous_symbol_mappings = collections.OrderedDict()
        self.ambiguous_symbol_mappings["?"] = "ACGTU-"
        self.ambiguous_symbol_mappings["N"] = "ACGTU"
        self.ambiguous_symbol_mappings["R"] = "AG"
        self.ambiguous_symbol_mappings["Y"] = "CTU"
        self.ambiguous_symbol_mappings["M"] = "AC"
        self.ambiguous_symbol_mappings["W"] = "ATU"
        self.ambiguous_symbol_mappings["S"] = "CG"
        self.ambiguous_symbol_mappings["K"] = "GTU"
        self.ambiguous_symbol_mappings["V"] = "ACG"
        self.ambiguous_symbol_mappings["H"] = "ACTU"
        self.ambiguous_symbol_mappings["D"] = "AGTU"
        self.ambiguous_symbol_mappings["B"] = "CGTU"

        self.polymorphic_symbol_mappings = collections.OrderedDict()

        # note reverse polarity here: from referenced to referencing
        self.additional_synonyms_map = collections.OrderedDict()
        self.additional_synonyms_map["N"] = "X"

        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())

        self.sa = dendropy.NUCLEOTIDE_STATE_ALPHABET
        self.num_total_states = len(self.expected_fundamental_state_symbols) + len(self.ambiguous_symbol_mappings) + len(self.polymorphic_symbol_mappings)

class ProteinStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.random_seed = random.randint(0, sys.maxsize)
        self.rng = random.Random(self.random_seed)

        self.expected_fundamental_state_symbols = [
                "A", "C", "D", "E", "F", "G", "H", "I",
                "K", "L", "M", "N", "P", "Q", "R", "S",
                "T",  "V", "W", "Y", "*", "-",
                ]
        self.ambiguous_symbol_mappings = collections.OrderedDict()
        self.ambiguous_symbol_mappings["?"] = "ACDEFGHIKLMNPQRSTVWY*-"
        self.ambiguous_symbol_mappings["B"] = "DN"
        self.ambiguous_symbol_mappings["Z"] = "EQ"
        self.ambiguous_symbol_mappings["X"] = "ACDEFGHIKLMNPQRSTVWY*"

        self.polymorphic_symbol_mappings = collections.OrderedDict()

        # note reverse polarity here: from referenced to referencing
        self.additional_synonyms_map = collections.OrderedDict()
        # self.additional_synonyms_map["N"] = "X"

        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())

        self.sa = dendropy.PROTEIN_STATE_ALPHABET
        self.num_total_states = len(self.expected_fundamental_state_symbols) + len(self.ambiguous_symbol_mappings) + len(self.polymorphic_symbol_mappings)

class BinaryStateAlphabetTest(
        StateAlphabetTester,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.random_seed = random.randint(0, sys.maxsize)
        self.rng = random.Random(self.random_seed)

        self.expected_fundamental_state_symbols = ["1", "0"]
        self.ambiguous_symbol_mappings = collections.OrderedDict()

        self.polymorphic_symbol_mappings = collections.OrderedDict()

        # note reverse polarity here: from referenced to referencing
        self.additional_synonyms_map = collections.OrderedDict()
        # self.additional_synonyms_map["N"] = "X"

        self.expected_polymorphic_state_symbols = list(self.polymorphic_symbol_mappings.keys())
        self.expected_ambiguous_state_symbols = list(self.ambiguous_symbol_mappings.keys())

        self.sa = dendropy.BINARY_STATE_ALPHABET
        self.num_total_states = len(self.expected_fundamental_state_symbols) + len(self.ambiguous_symbol_mappings) + len(self.polymorphic_symbol_mappings)


if __name__ == "__main__":
    unittest.main()
